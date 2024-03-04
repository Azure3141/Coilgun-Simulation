import numpy as np
import matplotlib.pyplot as plt
import run_step as run
import classes
import parameters

copper = parameters.copper

# Coils
# Material, inner radius, length, layers, wire diameter:
armature = classes.Coil(copper, 0.0125, 0.05, 1, 0.00129)

driver1 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)
driver2 = classes.Coil(copper, 0.015, 0.05, 2, 0.00129)
driver3 = classes.Coil(copper, 0.015, 0.05, 2, 0.00129)
driver4 = classes.Coil(copper, 0.015, 0.05, 2, 0.00129)
driver5 = classes.Coil(copper, 0.015, 0.05, 2, 0.00129)

# Stages
# Coil, position, trigger position, capacitance, charge voltage
stage1 = classes.Stage(driver1, 0.0, 0.0, 0.0024 * 2, 450)
stage2 = classes.Stage(driver2, 0.1, 0.1, 0.0024 * 2, 450)
stage3 = classes.Stage(driver3, 0.2, 0.2, 0.0024 * 2, 450)
stage4 = classes.Stage(driver4, 0.3, 0.3, 0.0024 * 2, 450)
stage5 = classes.Stage(driver5, 0.4, 0.4, 0.0024 * 2, 450)

# Projectile
# Coil, initial position, mass
projectile = classes.Projectile(armature, 0.01, 0.1)

stage_list = [stage1, stage2, stage3, stage4, stage5]
efficiency = []

runtime = 0.004
timestep = 0.0001

for step in np.arange(0, runtime, timestep):
    # Timestep, driver list, projectile
    output = run.run_turn(timestep, stage_list, projectile)

    driver_energy_sum = driver1.energy + driver2.energy + driver3.energy + driver4.energy + driver5.energy
    efficiency.append(armature.energy / (driver_energy_sum + 0.001))

time = np.arange(0, runtime, timestep)

# figp, p = plt.subplots()
# p.plot(time, projectile.position_list)
# p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
# figp.savefig("pos.png")

figv, v = plt.subplots()
v.plot(time, projectile.velocity_list)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figa, p = plt.subplots()
p.plot(time, projectile.acceleration_list)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("accel.png")

figi, i = plt.subplots()
i.plot(time, projectile.coil.currents_list)
for stage in range(len(stage_list)):
    i.plot(time, stage_list[stage].coil.currents_list)
i.set(xlabel='Time (s)', ylabel='Currents (A)', title='Currents')
i.legend(['Armature', 'Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5'])
figi.savefig("currents.png")

# figv, v = plt.subplots()
# v.plot(time, projectile.coil.voltage_list)
# for stage in range(len(stage_list)):
#     v.plot(time, stage_list[stage].coil.voltage_list)
# v.set(xlabel='Time (s)', ylabel='Voltages (V)', title='Voltages')
# v.legend(['Armature', 'Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5'])
# figv.savefig("voltages.png")

fign, n = plt.subplots()
n.plot(time, efficiency)
n.set(xlabel='Time (s)', ylabel='Efficiency', title='Coilgun Efficiency')
fign.savefig("efficiency.png")

fige, e = plt.subplots()
e.plot(time, projectile.coil.energy_list)
for stage in range(len(stage_list)):
    e.plot(time, stage_list[stage].coil.energy_list)
e.set(xlabel='Time (s)', ylabel='Energy (J)', title='Energy')
e.legend(['Armature', 'Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5'])
fige.savefig("energy.png")

figt, t = plt.subplots()
t.plot(time, projectile.coil.temperature_list)
for stage in range(len(stage_list)):
    t.plot(time, stage_list[stage].coil.temperature_list)
t.set(xlabel='Time (s)', ylabel='Temperature (K)', title='Temperature')
t.legend(['Armature', 'Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5'])
figt.savefig("temperature.png")


a1 = stage1.coil.r_inner - stage1.coil.wire_d / 2
a2 = stage1.coil.r_outer + stage1.coil.wire_d / 2
w = (a2 - a1) / parameters.stress_element_divisions

figs, (hoop, radial) = plt.subplots(nrows=2)
hoopc = hoop.contourf(time, np.arange(a1, a2, w), stage1.coil.hoop_stress_list.T, cmap='plasma')
hoopcb = figs.colorbar(hoopc, ax=hoop)
hoopcb.set_label("Stress (Pa)", loc='center')
hoop.set(xlabel='Time (s)', ylabel='Radial Distance', title='Hoop Stress (Pa)')

radialc = radial.contourf(time, np.arange(a1, a2, w), stage1.coil.radial_stress_list.T, cmap='plasma')
radialcb = figs.colorbar(radialc, ax=radial)
radialcb.set_label("Stress (Pa)", loc='center')
radial.set(xlabel='Time (s)', ylabel='Radial Distance', title='Radial Stress (Pa)')
figs.savefig("stress.png")



plt.show()
