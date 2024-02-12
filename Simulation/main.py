import math
import matplotlib.pyplot as plt
from numba import jit
import numpy as np
import scipy as sp
import inductance
from inductance import elliptics


mu0 = 4 * np.pi * pow(10, -7)

class Coil:
    def __init__(self, position, r_inner, length, layers, wire_d, current):
        self.pos = position # Refers to axial distance of center of coil

        self.r_inner = r_inner
        self.length = length
        self.layers = layers
        self.wire_d = wire_d
        self.current = current

        self.turns = round(length * layers / wire_d)
        self.r_outer = r_inner + (layers - 1) * wire_d
        self.r_avg = (r_inner + self.r_outer) / 2
        self.area = pow(r_inner, 2) * np.pi
        self.wire_length = self.turns * np.pi * self.r_avg * 2

        self.resistance = 1.68 * pow(10, -8) * self.wire_length / (math.pi * pow(self.wire_d / 2, 2))
        self.inductance = mu0 * self.area * pow(self.turns, 2) / length
        self.tau = self.inductance / self.resistance

        print("Turns", self.turns)
        print("Inductance", self.inductance, "H")
        print("Resistance", self.resistance, "Ohm")
        print("Time Constant", self.tau)


class Thing:
    def __init__(self, mass, position, velocity, coil):
        self.mass = mass
        self.pos = position # Refers to axial distance of center of coil
        self.velocity = velocity
        self.coil = coil


def mutual_inductance(coil1, coil2):
    c1 = coil1
    c2 = coil2

    M_grad = 0
    M = 0

    for i in range(int(c1.turns)):
        for j in range(int(c2.turns)):
            #pos_i = c1.pos + i * c1.wire_d - c1.length / 2
            #pos_j = c2.pos + j * c2.wire_d - c2.length / 2

            pos_i = c1.pos - c1.length / 2 + c1.wire_d / c1.layers * i
            pos_j = c2.pos - c1.length / 2 + c2.wire_d / c1.layers * j

            z = pos_i - pos_j

            k = 2*math.sqrt(coil1.r_avg * coil2.r_avg) / math.sqrt(pow((coil1.r_avg + coil2.r_avg), 2) + pow(z, 2))

            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            M_grad = M_grad + (mu0 * k * z / (4 * (1 - pow(k, 2)) * math.sqrt(coil1.r_avg * coil2.r_avg)) * (2 * (1 - pow(k, 2)) * k_k + (2 - pow(k, 2)) * k_e))
            M = M + (-mu0 * math.sqrt(coil1.r_avg * coil2.r_avg) * ((k - 2 / k) * k_k + 2 / k * k_e))

            #F1 = F1 + mu0 * k / (4*math.sqrt(coil1.r_avg * coil2.r_avg)) * ((2 - pow(k, 2)) / (1 - pow(k, 2)) * k_e - 2 * k_k)

    return M, M_grad


def run_turn(timestep_size, step):
    driver.pos = 0
    m_armature = projectile.mass
    driver.current = 1000 * (1 - pow(math.e, -step / driver.tau))
    I_driver = driver.current

    armature.pos = projectile.pos

    M_vec = mutual_inductance(driver, armature)
    M = M_vec[0]
    M_grad = M_vec[1]

    I_armature = -M / armature.inductance * I_driver
    a_armature = M_grad * I_driver * I_armature / m_armature

    # Increment positions and velocities
    v_armature = projectile.velocity + a_armature * timestep_size
    projectile.pos = armature.pos + v_armature * timestep_size

    projectile.velocity = v_armature
    armature.current = I_armature

    return projectile.pos, projectile.velocity, a_armature, driver.current, armature.current, M, M_grad, armature.pos


# Initial position, inner radius, length, layers, wire diameter, current
driver = Coil(0, 0.015, 0.063, 2, 0.0014, 7500)
armature = Coil(0.01, 0.0125, 0.06, 1, 0.0014, 0)

# Mass, position, velocity, coil
projectile = Thing(0.1, armature.pos, 0, armature)

runtime = 0.1
timestep_size = 0.0001

positions = []
velocities = []
accelerations = []
driver_current = []
armature_current = []
m_inductances = []
mg_inductances = []

a_positions = []

for i in np.arange(0, runtime, timestep_size):
    output = run_turn(timestep_size, i)
    positions.append(output[0])
    velocities.append(output[1])
    accelerations.append(output[2])
    driver_current.append(output[3])
    armature_current.append(output[4])
    m_inductances.append(output[5])
    mg_inductances.append(output[6])
    a_positions.append(output[7])




figp, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep_size), positions)
p.plot(np.arange(0, runtime, timestep_size), a_positions)
p.set(xlabel='time (s)', ylabel='position (m)', title='Projectile Position')
figp.savefig("pos.png")

figv, v = plt.subplots()
v.plot(np.arange(0, runtime, timestep_size), velocities)
v.set(xlabel='time (s)', ylabel='velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figp, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep_size), accelerations)
p.set(xlabel='time (s)', ylabel='acceleration (m/s^2)', title='Projectile Acceleration')
figp.savefig("accel.png")

# figmp, mp = plt.subplots()
# mp.plot(positions, m_inductances)
# mp.set(xlabel='position (m)', ylabel='mutual inductance (H)', title='Mutual Inductance vs Position')
# figmp.savefig("mp.png")
#
# figmgp, mgp = plt.subplots()
# mgp.plot(positions, mg_inductances)
# mgp.set(xlabel='position (m)', ylabel='mutual inductance gradient (H/m)', title='M_grad vs Position')
# figmgp.savefig("mgp.png")

figp, i = plt.subplots()
i.plot(np.arange(0, runtime, timestep_size), driver_current)
i.plot(np.arange(0, runtime, timestep_size), armature_current)
i.set(xlabel='time (s)', ylabel='current (A)', title='Driver and Armature Currents')
i.legend(['Driver', 'Armature'])
figp.savefig("current.png")

# figm, m = plt.subplots()
# m.plot(np.arange(0, runtime, timestep_size), m_inductances)
# m.set(xlabel='time (s)', ylabel='mutual inductance (H)', title='Mutual Inductance')
# figm.savefig("m_inductance.png")
#
# figmg, mg = plt.subplots()
# mg.plot(np.arange(0, runtime, timestep_size), mg_inductances)
# mg.set(xlabel='time (s)', ylabel='mutual inductance gradient (H/m)', title='Mutual Inductance Gradient')
# figmg.savefig("mg_inductance.png")




# mu = []
# mug = []
# distances = np.arange(-1, 1, 0.001)
# for i in distances:
#     armature.pos = i
#     M_vec = mutual_inductance(driver, armature)
#     M = M_vec[0]
#     M_grad = M_vec[1]
#
#     mu.append(M)
#     mug.append(M_grad)
#
#
# figp, mup = plt.subplots()
# mup.plot(distances, mu)
# mup.set(xlabel='displacement (m)', ylabel='mutual inductance (H)', title='Mutual Inductance')
# figp.savefig("mu.png")
#
# figp, mugp = plt.subplots()
# mugp.plot(distances, mug)
# mugp.set(xlabel='displacement (m)', ylabel='mutual inductance gradient (H/m)', title='Mutual Inductance Gradient')
# figp.savefig("mug.png")




plt.show()