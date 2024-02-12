import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# Constants
mu0 = 4 * np.pi * pow(10, -7)
ambient_temp = 300


class Material:
    def __init__(self, density, heat_capacity, resistivity, alpha):
        self.density = density
        self.heat_capacity = heat_capacity
        self.resistivity = resistivity
        self.alpha = alpha


class Coil:
    def __init__(self, type, material, position, r_inner, length, layers, wire_d, capacitance, initial_voltage, trigger_pos):
        self.type = type
        self.pos = position  # Refers to axial distance of coil center
        self.trigger_pos = trigger_pos
        self.vel = 0
        self.current = 0
        self.energy = 0
        self.temp = ambient_temp
        self.switch_step = 0
        self.switch = 0
        self.state = 0

        self.density = material.density
        self.heat_capacity = material.heat_capacity
        self.alpha = material.alpha
        self.resistivity = material.resistivity
        self.base_resistivity = self.resistivity

        self.mass = 0
        self.r_inner = r_inner
        self.length = length
        self.layers = layers
        self.wire_d = wire_d
        self.peak_current = 0
        self.initial_voltage = initial_voltage

        self.turns = round(length * layers / wire_d)
        self.r_outer = r_inner + (layers - 1) * wire_d
        self.r_avg = (r_inner + self.r_outer) / 2
        self.area = pow(r_inner, 2) * np.pi
        self.wire_length = self.turns * np.pi * self.r_avg * 2
        self.wire_area = math.pi * pow(self.wire_d / 2, 2)

        self.wire_mass = self.wire_area * self.wire_length * self.density
        self.thermal_mass = self.wire_mass * self.heat_capacity

        self.capacitance = capacitance
        self.inductance = mu0 * self.area * pow(self.turns, 2) / length
        self.resistance = self.resistivity * self.wire_length / self.wire_area
        self.wire_resistance = self.resistance

        if(self.type == "Driver"):
            self.resistance = 2 * math.sqrt(self.inductance / self.capacitance)

        self.added_resistance = self.resistance - self.wire_resistance
        self.tau = self.inductance / self.resistance

        print("Wire Length", self.wire_length)
        print("Turns", self.turns)
        print("Inductance", self.inductance, "H")
        #print("Required Resistance", self.resistance, "Ohm")
        print("Added Resistance", self.added_resistance, "Ohm")
        #print("Time Constant", self.tau)


def mutual_inductance(coil1, coil2):
    M = 0
    for i in range(int(coil1.turns)):
        for j in range(int(coil2.turns)):
            pos_i = coil1.pos - coil1.length / 2 + coil1.wire_d / coil1.layers * i
            pos_j = coil2.pos - coil1.length / 2 + coil2.wire_d / coil1.layers * j
            z = pos_i - pos_j

            k = 2*math.sqrt(coil1.r_avg * coil2.r_avg) / math.sqrt(pow((coil1.r_avg + coil2.r_avg), 2) + pow(z, 2))
            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            M = M + (-mu0 * math.sqrt(coil1.r_avg * coil2.r_avg) * ((k - 2 / k) * k_k + 2 / k * k_e))

    return M


def run_turn(timestep, step, stages, armature):

    F = 0
    # Run switching
    new_step = step - armature.switch_step
    switch = armature.switch


    for stage in stages:
        driver = stage
        M = mutual_inductance(driver, armature)

        if(armature.pos >= driver.trigger_pos and armature.switch < len(stages)):
            armature.switch_step = step
            armature.switch += 1
            driver.state = 1
            print("Switched stage", switch + 1)

        # Run currents
        #I_driver = coil1.peak_current * (1 - pow(math.e, -new_step / coil1.tau))
        # I_deriv = coil1.peak_current * pow(math.e, -timestep / coil1.tau) / coil1.tau

        I_driver = driver.initial_voltage / driver.inductance * new_step * \
                   pow(math.e, -driver.resistance / (2 * driver.inductance) * new_step) * driver.state

        driver.current = I_driver
        I_armature = -M / armature.inductance * driver.current * (1 - pow(math.e, -new_step / armature.tau))
        armature.current = I_armature

        # Run thermals
        driver_power = pow(I_driver, 2) * driver.resistance
        armature_power = pow(I_armature, 2) * armature.resistance
        driver.temp = driver.temp + driver_power / driver.thermal_mass * timestep
        armature.temp = armature.temp + armature_power / armature.thermal_mass * timestep

        # Run resistivity
        # coil1.resistivity = coil1.base_resistivity + coil1.alpha * (coil1.temp - ambient_temp)
        # coil1.resistance = coil1.resistivity * coil1.wire_length / coil1.wire_area
        armature.resistivity = armature.base_resistivity + armature.alpha * (armature.temp - ambient_temp)
        armature.resistance = armature.resistivity * armature.wire_length / armature.wire_area
        armature.tau = armature.inductance / armature.resistance

        # Energy balance
        driver.energy = driver.energy + driver_power * timestep
        armature.energy = pow(armature.vel, 2) * 0.5 * armature.mass

        # Solve magnetic force
        K = 20
        m = K
        a = driver.length
        b = armature.length
        r1 = driver.r_avg
        r2 = armature.r_avg
        F_sum = 0

        for i in range(2*K+1):
            for j in range(2*m+1):
                g = -K + i
                l = -m + j
                c = driver.pos - armature.pos
                z = c - a * g / (2 * K + 1) + b * l / (2 * m + 1)

                k = math.sqrt(4 * r1 * r2 / (pow(r1 + r2, 2) + pow(z, 2)))
                k_k = sp.special.ellipk(pow(k, 2))
                k_e = sp.special.ellipe(pow(k, 2))

                F_sum = F_sum + (mu0 * z * I_armature * I_driver * k / (4 * math.sqrt(r1 * r2))
                                 * ((2 - pow(k, 2)) / (1 - pow(k, 2)) * k_e - 2*k_k))

        F_coeff = driver.turns * armature.turns / ((2 * K + 1) * (2 * m + 1))
        F += F_coeff * F_sum

    # Run kinematics
    a_armature = F / armature.mass
    armature.vel = armature.vel + a_armature * timestep
    armature.pos = armature.pos + armature.vel * timestep

    return armature.pos, armature.vel, a_armature, driver.current, armature.current, M, driver.temp, armature.temp, driver.energy, armature.energy


# Density, heat capacity, resistivity, alpha
copper = Material(8960, 385, 1.68 * pow(10, -8), 0.00382 * pow(10, -8))

# Type, initial position, material, inner radius, length, layers, wire diameter, capacitance, voltage, trigger position
# driver = Coil("Driver", copper, 0, 0.015, 0.05, 3, 0.001024, 0.0048)
stage1 = Coil("Driver", copper, 0.1, 0.015, 0.05, 1, 0.001024, 0.0048, 450, 0)
stage2 = Coil("Driver", copper, 0.2, 0.015, 0.05, 2, 0.001024, 0.0048, 450, 0.1)
stage3 = Coil("Driver", copper, 0.3, 0.015, 0.05, 3, 0.001024, 0.0048, 450, 0.2)
stage4 = Coil("Driver", copper, 0.4, 0.015, 0.05, 4, 0.001024, 0.0048, 450, 0.3)
stage5 = Coil("Driver", copper, 0.5, 0.015, 0.05, 5, 0.001024, 0.0048, 450, 0.4)

stages = [stage1, stage2, stage3, stage4, stage5]
armature = Coil("Proj", copper, 0.01, 0.0125, 0.05, 1, 0.001024, 0, 0, 0)
armature.mass = 0.1



# Stage positions and trigger positions
positions = [0, 0, 0, 0, 0]
triggers = [0, 0.1, 0.2, 0.3, 0.4]

# Capacitance, charge voltage
capacitors = [0.0048, 450]

# Inner radius, length, layers, wire diameter
geometries = [0.015, 0.05, 3, 0.001024]


#driver_positions = [0, 0.1, 0.2, 0.3, 0.4]
#trigger_list = [0, 0.1, 0.2, 0.3, 0.4, 100]



runtime = 0.05
timestep = 0.0001

positions = []
velocities = []
accelerations = []
driver_current = []
armature_current = []
m_inductances = []
driver_temp = []
armature_temp = []
driver_energy = []
armature_energy = []
efficiency = []

for i in np.arange(0, runtime, timestep):
    # Timestep, step, coil1, coil2, driver positions, trigger positions
    output = run_turn(timestep, i, stages, armature)
    positions.append(output[0])
    velocities.append(output[1])
    accelerations.append(output[2])
    driver_current.append(output[3])
    armature_current.append(output[4])
    m_inductances.append(output[5])
    driver_temp.append(output[6])
    armature_temp.append(output[7])
    driver_energy.append(output[8])
    armature_energy.append(output[9])
    efficiency.append(output[9] / (output[8] + 0.001))


figp, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep), positions)
p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
figp.savefig("pos.png")

figv, v = plt.subplots()
v.plot(np.arange(0, runtime, timestep), velocities)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figa, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep), accelerations)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("accel.png")

figi, i = plt.subplots()
i.plot(np.arange(0, runtime, timestep), driver_current)
i.plot(np.arange(0, runtime, timestep), armature_current)
i.set(xlabel='Time (s)', ylabel='Current (A)', title='Driver and Armature Currents')
i.legend(['Driver', 'Armature'])
figi.savefig("current.png")

figth, th = plt.subplots()
th.plot(np.arange(0, runtime, timestep), driver_temp)
th.plot(np.arange(0, runtime, timestep), armature_temp)
th.set(xlabel='Time (s)', ylabel='Temperature (K)', title='Driver and Armature Temperatures')
th.legend(['Driver', 'Armature'])
figth.savefig("temperatures.png")

fige, e = plt.subplots()
e.plot(np.arange(0, runtime, timestep), driver_energy)
e.plot(np.arange(0, runtime, timestep), armature_energy)
e.set(xlabel='Time (s)', ylabel='Energy (J)', title='Energy Balance')
e.legend(['Coil Electrical Energy', 'Projectile Kinetic Energy'])
fige.savefig("energy.png")

figeff, eff = plt.subplots()
eff.plot(np.arange(0, runtime, timestep), efficiency)
eff.set(xlabel='Time (s)', ylabel='Efficiency', title='Energy Efficiency')
figeff.savefig("efficiency.png")

figm, m = plt.subplots()
m.plot(np.arange(0, runtime, timestep), m_inductances)
m.set(xlabel='Time (s)', ylabel='Mutual Inductance (H)', title='Mutual Inductance')
figm.savefig("m_inductance.png")


# figmg, mg = plt.subplots()
# mg.plot(np.arange(0, time, timestep_size), mg_inductances)
# mg.set(xlabel='time (s)', ylabel='mutual inductance gradient (H/m)', title='Mutual Inductance Gradient')
# figmg.savefig("mg_inductance.png")
#
# figmgp, mgp = plt.subplots()
# mgp.plot(positions, mg_inductances)
# mgp.set(xlabel='position (m)', ylabel='mutual inductance gradient (H/m)', title='M_grad vs Position')
# figmgp.savefig("mgp.png")


# mu = []
# mug = []
# distances = np.arange(-1, 1, 0.01)
# for i in distances:
#     armature.pos = i
#     M_vec = mutual_inductance(driver, armature)
#     M = M_vec[0]
#     M_grad = M_vec[1]
#
#     mu.append(M)
#     mug.append(M_grad)


# figp, mup = plt.subplots()
# mup.plot(distances, mu)
# mup.set(xlabel='Displacement (m)', ylabel='mutual inductance (H)', title='Mutual Inductance')
# figp.savefig("mu.png")
#
# figp, mugp = plt.subplots()
# mugp.plot(distances, mug)
# mugp.set(xlabel='displacement (x)', ylabel='mutual inductance gradient (H/m)', title='Mutual Inductance Gradient')
# figmugp.savefig("mug.png")


plt.show()