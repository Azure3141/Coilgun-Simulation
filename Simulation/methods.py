import math
import scipy as sp
import parameters

def mutual_inductance(stage, projectile):
    coil1 = stage.coil
    coil2 = projectile.coil
    M = 0
    for i in range(int(coil1.turns)):
        for j in range(int(coil2.turns)):
            pos_i = stage.pos - coil1.length / 2 + coil1.wire_d / coil1.layers * i
            pos_j = projectile.pos - coil1.length / 2 + coil2.wire_d / coil1.layers * j
            z = pos_i - pos_j

            k = 2*math.sqrt(coil1.r_avg * coil2.r_avg) / math.sqrt(pow((coil1.r_avg + coil2.r_avg), 2) + pow(z, 2))
            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            M = M + (-parameters.mu0 * math.sqrt(coil1.r_avg * coil2.r_avg) * ((k - 2 / k) * k_k + 2 / k * k_e))

    return M


def solve_currents(stage, projectile, driver_step, armature_step):
    armature = projectile.coil
    driver = stage.coil
    M = mutual_inductance(stage, projectile)
    driver.current = stage.voltage / driver.inductance * driver_step * \
                          pow(math.e, -driver.resistance / (2 * driver.inductance) * driver_step)
    armature.current = -M / armature.inductance * driver.current * (1 - pow(math.e, -armature_step / armature.tau))


def increment_time(timestep, coil):
    coil.time_on += timestep


def solve_energy(coil, timestep):
    coil_power = pow(coil.current, 2) * coil.resistance
    coil.temp = coil.temp + coil_power / coil.thermal_mass * timestep
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


