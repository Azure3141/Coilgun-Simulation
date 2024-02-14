import numpy as np
import math

import run_step
import material
import coil

# Constants
mu0 = 4 * np.pi * pow(10, -7)
ambient_temp = 300

# Density, heat capacity, resistivity, alpha
copper = material.Material(8960, 385, 1.68 * pow(10, -8), 0.00382 * pow(10, -8))

# Type, initial position, material, inner radius, length, layers, wire diameter, capacitance
driver = coil.Coil("Driver", copper, 0, 0.015, 0.05, 3, 0.001024, 0.0048)
armature = coil.Coil("Proj", copper, 0.01, 0.0125, 0.05, 1, 0.001024, 0)

driver_positions = [0, 0.1, 0.2, 0.3, 0.4]
trigger_list = [0, 0.1, 0.2, 0.3, 0.4, 100]

armature.mass = 0.1
driver.initial_voltage = 450

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
    output = run_step.run_turn(timestep, i, driver, armature, trigger_list)
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