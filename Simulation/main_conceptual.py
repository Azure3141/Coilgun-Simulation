import numpy as np
import matplotlib.pyplot as plt
import run_step as run
import classes
import parameters

copper = parameters.copper

# # Coils
# # Material, inner radius, length, layers, wire diameter:
# armature = classes.Coil(copper, 0.33 / 2, 3, 1, 0.025)
#
# # Projectile
# # Coil, initial position, mass
# projectile = classes.Projectile(armature, 0.1, 1250)
#
# stage_count = 20
# stage_list = []
#
# # Material, inner radius, length, layers, wire diameter:
# drivers = [classes.Coil(copper, 0.38 / 2, 3, 3, 0.025) for i in range(stage_count)]
#
# # Coil, position, trigger position, capacitance, charge voltage
# stages = [classes.Stage(drivers[j], 3.0 * j, 3.0 * j, 0.01, 50000) for j in range(stage_count)]
#
# for i in range(stage_count):
#     stage_list.append(stages[i])



stage_count = 20
stage_list = []

# Material, inner radius, length, layers, wire diameter:
armature = classes.Coil(copper, 0.33, 2, 1, 0.025)
# Coil, initial position, mass
projectile = classes.Projectile(armature, 0.05, 1250)

# Material, inner radius, length, layers, wire diameter:
drivers = [classes.Coil(copper, 0.35, 2, 1, 0.025) for i in range(stage_count)]
# Coil, position, trigger position, capacitance, charge voltage
stages = [classes.Stage(drivers[j], 2.0 * j, 2.0 * j, 0.02, 50000) for j in range(stage_count)]

stages[0].capacitance = 0.075
stages[0].layers = 3
stages[1].capacitance = 0.05
stages[1].layers = 2

for i in range(stage_count):
    stage_list.append(stages[i])


efficiency = []
driver_energy_sum = 0.
driver_total_energy = []

runtime = 0.25
timestep = 0.0001

for step in np.arange(0, runtime, timestep):
    # Timestep, driver list, projectile
    output = run.run_turn(timestep, stage_list, projectile)

    driver_energy_sum = sum([drivers[k].energy for k in range(stage_count)])
    driver_total_energy.append(driver_energy_sum)
    efficiency.append(armature.energy / (driver_energy_sum + 0.001))

time = np.arange(0, runtime, timestep)

figp, p = plt.subplots()
p.plot(time, projectile.position_list)
p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
figp.savefig("conceptual_pos.png")

figv, v = plt.subplots()
v.plot(time, projectile.velocity_list)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("conceptual_vel.png")

figa, p = plt.subplots()
p.plot(time, projectile.acceleration_list)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("conceptual_accel.png")

figi, i = plt.subplots()
i.plot(time, projectile.coil.currents_list)
ilegend = ["Armature"]
for stage in range(len(stage_list)):
    i.plot(time, stage_list[stage].coil.currents_list)
    ilegend.append(("Driver " + str(stage + 1)))
i.set(xlabel='Time (s)', ylabel='Currents (A)', title='Currents')
i.legend(ilegend)
figi.savefig("conceptual_currents.png")

# figv, v = plt.subplots()
# v.plot(time, projectile.coil.voltage_list)
# for stage in range(len(stage_list)):
#     v.plot(time, stage_list[stage].coil.voltage_list)
# v.set(xlabel='Time (s)', ylabel='Voltages (V)', title='Voltages')
# v.legend(['Armature', 'Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5'])
# figv.savefig("conceptual_voltages.png")

fign, n = plt.subplots()
n.plot(time, efficiency)
n.set(xlabel='Time (s)', ylabel='Efficiency', title='Coilgun Efficiency')
fign.savefig("conceptual_efficiency.png")

fige, e = plt.subplots()
e.plot(time, projectile.coil.energy_list)
e.plot(time, driver_total_energy)
elegend = ["Armature", "Total Energy"]
for stage in range(len(stage_list)):
    e.plot(time, stage_list[stage].coil.energy_list)
    elegend.append(("Driver " + str(stage + 1)))
e.set(xlabel='Time (s)', ylabel='Energy (J)', title='Energy')
e.legend(elegend)
fige.savefig("conceptual_energy.png")

figt, t = plt.subplots()
t.plot(time, projectile.coil.temperature_list)
tlegend = ["Armature"]
for stage in range(len(stage_list)):
    t.plot(time, stage_list[stage].coil.temperature_list)
    tlegend.append(("Driver " + str(stage + 1)))
t.set(xlabel='Time (s)', ylabel='Temperature (K)', title='Temperature')
t.legend(tlegend)
figt.savefig("conceptual_temperature.png")


# a1 = stages[0].coil.r_inner - stages[0].coil.wire_d / 2
# a2 = stages[0].coil.r_outer + stages[0].coil.wire_d / 2
# w = (a2 - a1) / parameters.stress_element_divisions
#
# figs, (hoop, radial) = plt.subplots(nrows=2)
# hoopc = hoop.contourf(time, np.arange(a1, a2, w), stages[0].coil.hoop_stress_list.T, cmap='plasma')
# hoopcb = figs.colorbar(hoopc, ax=hoop)
# hoopcb.set_label("Stress (Pa)", loc='center')
# hoop.set(xlabel='Time (s)', ylabel='Radial Distance', title='Hoop Stress (Pa)')
#
# radialc = radial.contourf(time, np.arange(a1, a2, w), stages[0].coil.radial_stress_list.T, cmap='plasma')
# radialcb = figs.colorbar(radialc, ax=radial)
# radialcb.set_label("Stress (Pa)", loc='center')
# radial.set(xlabel='Time (s)', ylabel='Radial Distance', title='Radial Stress (Pa)')
# figs.savefig("conceptual_stress.png")



plt.show()
