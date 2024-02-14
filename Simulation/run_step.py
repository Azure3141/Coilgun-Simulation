import math
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

import magnetic_force
import inputs


def run_turn(timestep, step, drivers, armature):
    e


def run_turn2(timestep, step, drivers, armature):
    coil1 = drivers
    coil2 = armature
    # Run switching
    new_step = step - coil1.switch_step
    switch = min(coil1.switch, len(trigger_positions) - 1)
    M = mutual_inductance(coil1, coil2)
    if(coil2.pos >= trigger_positions[switch]):
        coil1.switch_step = step
        coil1.switch += 1
        coil1.pos = inputs.driver_positions[switch]
        coil1.temp = inputs.ambient_temp
        print("Switched stage", switch + 1)

    # Run currents
    #I_driver = coil1.peak_current * (1 - pow(math.e, -new_step / coil1.tau))
    # I_deriv = coil1.peak_current * pow(math.e, -timestep / coil1.tau) / coil1.tau

    I_driver = coil1.initial_voltage / coil1.inductance * new_step * pow(math.e, -coil1.resistance / (2 * coil1.inductance) * new_step)

    coil1.current = I_driver
    I_armature = -M / coil2.inductance * coil1.current * (1 - pow(math.e, -new_step / coil2.tau))
    coil2.current = I_armature

    # Run thermals
    driver_power = pow(I_driver, 2) * coil1.resistance
    armature_power = pow(I_armature, 2) * coil2.resistance
    coil1.temp = coil1.temp + driver_power / coil1.thermal_mass * timestep
    coil2.temp = coil2.temp + armature_power / coil2.thermal_mass * timestep

    # Run resistivity
    # coil1.resistivity = coil1.base_resistivity + coil1.alpha * (coil1.temp - ambient_temp)
    # coil1.resistance = coil1.resistivity * coil1.wire_length / coil1.wire_area
    coil2.resistivity = coil2.base_resistivity + coil2.alpha * (coil2.temp - inputs.ambient_temp)
    coil2.resistance = coil2.resistivity * coil2.wire_length / coil2.wire_area
    coil2.tau = coil2.inductance / coil2.resistance

    # Energy balance
    driver.energy = driver.energy + driver_power * timestep
    armature.energy = pow(armature.vel, 2) * 0.5 * armature.mass

    F = magnetic_force(driver, armature)

    # Run kinematics
    a_armature = F / armature.mass
    armature.vel = armature.vel + a_armature * timestep
    armature.pos = armature.pos + armature.vel * timestep

    return armature.pos, armature.vel, a_armature, driver.current, armature.current, M, driver.temp, armature.temp, driver.energy, armature.energy