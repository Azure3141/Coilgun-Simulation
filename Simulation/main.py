import numpy as np
import matplotlib.pyplot as plt
import run_step as run
import classes
import parameters

copper = parameters.copper

# Coils
# Material, inner radius, length, layers, wire diameter):
armature = classes.Coil(copper, 0.0125, 0.05, 1, 0.00129)

driver1 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)
driver2 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)
driver3 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)
driver4 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)
driver5 = classes.Coil(copper, 0.015, 0.05, 3, 0.00129)

# Stages
# Coil, position, trigger position, capacitance, voltage
stage1 = classes.Stage(driver1, 0.0, 0.0, 0.0048, 450)
stage2 = classes.Stage(driver2, 0.1, 0.1, 0.0048, 450)
stage3 = classes.Stage(driver3, 0.2, 0.2, 0.0048, 450)
stage4 = classes.Stage(driver4, 0.3, 0.3, 0.0048, 450)
stage5 = classes.Stage(driver5, 0.4, 0.4, 0.0048, 450)

# Projectile
# Coil, position, mass
projectile = classes.Projectile(armature, 0.01, 0.1)


runtime = 0.02
timestep = 0.0001

stage_list = [stage1, stage2, stage3, stage4, stage5]
velocities = []
accelerations = []
positions = []
armature_currents = []

armature_energy = []
driver_energy = []
efficiency = []

for step in np.arange(0, runtime, timestep):
    # Timestep, step, driver list, projectile
    output = run.run_turn(timestep, stage_list, projectile)

    velocities.append(output[0].vel)
    accelerations.append(output[0].accel)
    positions.append(output[0].pos)
    armature_currents.append(output[0].coil.current)

    armature_energy.append(armature.energy)
    driver_energy.append(driver1.energy + driver2.energy + driver3.energy + driver4.energy + driver5.energy)
    efficiency.append(armature.energy / (driver1.energy + driver2.energy + driver3.energy + driver4.energy + driver5.energy))


time = np.arange(0, runtime, timestep)

figp, p = plt.subplots()
p.plot(time, positions)
p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
figp.savefig("pos.png")

figv, v = plt.subplots()
v.plot(time, velocities)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figa, p = plt.subplots()
p.plot(time, accelerations)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("accel.png")

figi, i = plt.subplots()
i.plot(time, driver1.currents_list, time, driver2.currents_list, time, driver3.currents_list, time, driver4.currents_list, time, driver5.currents_list, time, armature_currents)
i.set(xlabel='Time (s)', ylabel='Currents (A)', title='Currents')
i.legend(['Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5', 'Armature'])
figi.savefig("currents.png")

fign, n = plt.subplots()
n.plot(time, efficiency)
n.set(xlabel='Time (s)', ylabel='Efficiency', title='Coilgun Efficiency')
fign.savefig("efficiency.png")

fige, e = plt.subplots()
e.plot(time, driver1.energy_list, time, driver2.energy_list, time, driver3.energy_list, time, driver4.energy_list, time, driver5.energy_list, time, armature.energy_list)
e.set(xlabel='Time (s)', ylabel='Energy (J)', title='Energy')
e.legend(['Driver 1', 'Driver 2', 'Driver 3', 'Driver 4', 'Driver 5', 'Armature'])
fige.savefig("Energy.png")


plt.show()
