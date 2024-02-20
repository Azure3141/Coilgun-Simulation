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
payload = classes.Projectile(armature, 0.01, 0.1)


runtime = 0.01
timestep = 0.0001

stage_list = [stage1, stage2, stage3, stage4, stage5]
velocities = []
accelerations = []
positions = []

for step in np.arange(0, runtime, timestep):
    # Timestep, step, driver list, projectile
    output = run.run_turn(timestep, stage_list, payload)
    velocities.append(output.vel)
    accelerations.append(output.accel)
    positions.append(output.pos)


figp, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep), positions)
p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
figp.savefig("vel.png")

figv, v = plt.subplots()
v.plot(np.arange(0, runtime, timestep), velocities)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figa, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep), accelerations)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("accel.png")

plt.show()

