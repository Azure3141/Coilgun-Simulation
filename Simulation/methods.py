import math
import numpy as np
import scipy as sp
import parameters

def mutual_inductance(stage, projectile):
    coil1 = stage.coil
    coil2 = projectile.coil
    M = 0

    for l1 in range(coil1.layers):
        r1_eff = coil1.r_inner + l1 * coil1.wire_d
        for l2 in range(coil2.layers):
            r2_eff = coil2.r_inner + l2 * coil2.wire_d
            for i in range(coil1.layer_turns[l1]):
                for j in range(coil2.layer_turns[l2]):
                    pos_i = stage.pos - coil1.length / 2 + coil1.wire_d * i
                    pos_j = projectile.pos - coil2.length / 2 + coil2.wire_d * j
                    z = pos_i - pos_j

                    k = 2 * math.sqrt(r1_eff * r2_eff) / math.sqrt(pow((r1_eff + r2_eff), 2) + pow(z, 2))
                    k_k = sp.special.ellipk(pow(k, 2))
                    k_e = sp.special.ellipe(pow(k, 2))

                    M += -parameters.mu0 * math.sqrt(r1_eff * r2_eff) * ((k - 2 / k) * k_k + 2 / k * k_e)
    return M


def solve_driver_current(stage, driver_step):
    driver = stage.coil

    # driver.current = stage.initial_voltage / driver.inductance * driver_step * \
    #                       pow(math.e, -driver.resistance / (2 * driver.inductance) * driver_step)

    # Natural response frequency
    omega = math.sqrt(1/(driver.inductance * stage.capacitance) - pow(driver.resistance, 2) / (4 * pow(driver.inductance, 2)))

    # Diode becomes active when capacitor voltage reaches 0
    if(driver.voltage < 0):
        # RL current model
        driver.current = driver.peak_current * pow(math.e, -1 / driver.tau * (driver_step - driver.diode_step))
        driver.voltage = -0.1
    else:
        # RLC current model
        driver.current = (stage.capacitance * stage.initial_voltage * pow(math.e, -driver_step/(2 * driver.tau)) *
                          (4 * pow(omega, 2) + pow(1 / driver.tau, 2)) * math.sin(driver_step * omega)) / (4 * omega)

        # Capacitor voltage
        driver.voltage = stage.initial_voltage * pow(math.e, -driver_step/(2 * driver.tau)) * \
                         (1 / (2 * driver.tau * omega) * math.sin(omega * driver_step) + math.cos(omega * driver_step))

        driver.peak_current = driver.current
        driver.diode_step = driver_step

def solve_armature_current(stage, projectile, armature_step):
    armature = projectile.coil
    driver = stage.coil
    M = mutual_inductance(stage, projectile)

    # transfer_eff = pow(M, 2) / (armature.inductance * driver.inductance)
    transfer_eff = 1
    armature_current = transfer_eff * -M / armature.inductance * driver.current * (1 - pow(math.e, -armature_step / armature.tau))
    return armature_current

def solve_stresses(coil):
    nu = coil.material.poisson
    a1 = coil.r_inner - coil.wire_d / 2
    a2 = coil.r_outer + coil.wire_d / 2
    alpha = a2 / a1

    w = (a2 - a1) / parameters.stress_element_divisions
    hoop_stresses = []
    radial_stresses = []

    for r in np.arange(a1, a2, w):
        rho = r / a1
        J = coil.current / coil.wire_area
        B1 = parameters.mu0 * J * a1 * (alpha - 1)
        B2 = 0

        K = (alpha * B1 - B2) * J * a1 / (alpha - 1)
        M = (B2 - B1) * J * a1 / (alpha - 1)

        sigma_theta = K * (2 + nu) / (3 * (alpha + 1)) * (pow(alpha, 2) + alpha + 1 + pow(alpha, 2) / pow(rho, 2) - rho * (1 + 2 * nu) * (alpha + 1) / (2 + nu)) - M * (3 + nu) / 8 * (pow(alpha, 2) + 1 + pow(alpha, 2) / pow(rho, 2) - (1 + 3 * nu) / (3 + nu) * pow(rho, 2))
        sigma_r = K * (2 + nu) / (3 * (alpha + 1)) * (pow(alpha, 2) + alpha + 1 - pow(alpha, 2) / pow(rho, 2) - (alpha + 1) * rho) - M * (3 + nu) / 8 * (pow(alpha, 2) + 1 - pow(alpha, 2) / pow(rho, 2) - pow(rho, 2))

        hoop_stresses.append(sigma_theta)
        radial_stresses.append(sigma_r)

    coil.hoop_stress_list = np.vstack((coil.hoop_stress_list, hoop_stresses))
    coil.radial_stress_list = np.vstack((coil.radial_stress_list, radial_stresses))


def increment_time(timestep, coil):
    coil.time_on += timestep


def solve_energy(coil, timestep):
    coil_power = pow(coil.current, 2) * coil.resistance
    coil.energy = coil.energy + coil_power * timestep


def solve_thermals(coil, timestep):
    coil_power = pow(coil.current, 2) * coil.resistance
    coil.temp = coil.temp + coil_power / coil.thermal_mass * timestep


def solve_resistance(coil):
    coil.resistivity = coil.material.resistivity + coil.material.alpha * (coil.temp - parameters.ambient_temp)
    coil.resistance = coil.resistivity * coil.wire_length / coil.wire_area
    coil.tau = coil.inductance / coil.resistance


def magnetic_force(stage, projectile):
    coil1 = stage.coil
    coil2 = projectile.coil
    K = 20
    m = K
    a = coil1.length
    b = coil2.length
    F_sum = 0

    for i in range(2 * K + 1):
        for j in range(2 * m + 1):
            g = -K + i
            l = -m + j
            c = stage.pos - projectile.pos
            z = c - a * g / (2 * K + 1) + b * l / (2 * m + 1)

            k = math.sqrt(4 * coil1.r_avg * coil2.r_avg / (pow(coil1.r_avg + coil2.r_avg, 2) + pow(z, 2)))
            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            F_sum = F_sum + (parameters.mu0 * z * coil1.current * coil2.current * k / (4 * math.sqrt(coil1.r_avg * coil2.r_avg))
                             * ((2 - pow(k, 2)) / (1 - pow(k, 2)) * k_e - 2 * k_k))

    F_coeff = coil1.turns * coil2.turns / ((2 * K + 1) * (2 * m + 1))
    F = F_coeff * F_sum
    return F


def solve_kinematics(projectile, F, timestep):
    armature = projectile.coil

    projectile.accel = F / projectile.mass
    projectile.vel = projectile.vel + projectile.accel * timestep
    projectile.pos = projectile.pos + projectile.vel * timestep
    armature.pos = projectile.pos

    projectile.kinetic_energy = pow(projectile.vel, 2) * 0.5 * projectile.mass
    armature.energy = projectile.kinetic_energy
    projectile.energy_list.append(projectile.kinetic_energy)
    projectile.acceleration_list.append(projectile.accel)
    projectile.velocity_list.append(projectile.vel)
    projectile.position_list.append(projectile.pos)
