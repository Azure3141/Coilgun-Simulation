import numpy as np
import matplotlib.pyplot as plt
import run_step as run
import classes
import parameters

from scipy import optimize

copper = parameters.copper



def runSimulation(inputs):

    length, layers, wire_d = inputs

    layers = round(layers)

    # Coils
    # Material, inner radius, length, layers, wire diameter):
    driver1 = classes.Coil(copper, 0.015, length, layers, wire_d)
    armature = classes.Coil(copper, 0.0125, 0.1, 1, 0.001024)

    # Stages
    # Coil, position, trigger position, capacitance, voltage
    stage1 = classes.Stage(driver1, 0.0, 0.0, 0.003 * 2, 1000)

    # Projectile
    # Coil, position, mass
    projectile = classes.Projectile(armature, 0.01, 0.1)


    stage_list = [stage1]

    runtime = 0.01
    timestep = 0.0001

    for step in np.arange(0, runtime, timestep):
        # Timestep, step, driver list, projectile
        output = run.run_turn(timestep, stage_list, projectile)

    max_velocity = max(projectile.velocity_list)
    print(max_velocity, "m/s")

    del projectile
    del armature
    del driver1
    del stage1

    return 1 / max_velocity


guess = [0.05, 2, 0.00129]
constraints = [[0, 0.5], [1, 10], [0, 0.1]]
#result = optimize.minimize(runSimulation, guess, method="SLSQP" ,bounds=constraints)
result = optimize.differential_evolution(runSimulation, bounds=constraints)

print(result.x)

# runSimulation(guess)

# runtime = 0.01
# timestep = 0.0001
# time = np.arange(0, runtime, timestep)
#
# figv, v = plt.subplots()
# v.plot(time, projectile.velocity_list)
# v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
# figv.savefig("vel.png")
#
# plt.show()
