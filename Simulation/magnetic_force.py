import math
import scipy as sp
import inputs


def mutual_inductance(coil1, coil2):
    M = 0
    for i in range(int(coil1.turns)):
        for j in range(int(coil2.turns)):
            pos_i = coil1.pos - coil1.length / 2 + coil1.wire_d / coil1.layers * i
            pos_j = coil2.pos - coil1.length / 2 + coil2.wire_d / coil1.layers * j
            z = pos_i - pos_j

            k = 2*math.sqrt(coil1.r_avg * coil2.r_avg) / math.sqrt(pow((coil1.r_avg + coil2.r_avg), 2) + pow(z, 2))
            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            M = M + (-inputs.mu0 * math.sqrt(coil1.r_avg * coil2.r_avg) * ((k - 2 / k) * k_k + 2 / k * k_e))

    return M


def run_current(step, driver, armature):
    M = mutual_inductance(driver, armature)
    driver.current = driver.coil.initial_voltage / driver.inductance * step * pow(math.e, -driver.resistance / (2 * driver.inductance) * step)
    armature.current = -M / armature.inductance * driver.current * (1 - pow(math.e, -step / armature.tau))


def run_thermals(coil, timestep):
    coil_power = pow(coil.current, 2) * coil.resistance
    coil.temp = coil.temp + coil_power / coil.thermal_mass * timestep


def run_energy(coil, type, timestep):
    coil_power = pow(coil.current, 2) * coil.resistance
    coil.temp = coil.temp + coil_power / coil.thermal_mass * timestep
    
    if(type == 'kinetic'):
        coil.energy = pow(coil.vel, 2) * 0.5 * coil.mass
    elif(type == 'thermal'):
        coil.energy = coil.energy + coil_power * timestep


def run_resistance(coil):
    coil.resistivity = coil.base_resistivity + coil.alpha * (coil.temp - inputs.ambient_temp)
    coil.resistance = coil.resistivity * coil.wire_length / coil.wire_area
    coil.tau = coil.inductance / coil.resistance


def magnetic_force(coil1, coil2):
    K = 20
    m = K
    a = coil1.length
    b = coil2.length
    F_sum = 0

    for i in range(2 * K + 1):
        for j in range(2 * m + 1):
            g = -K + i
            l = -m + j
            c = coil1.pos - coil2.pos
            z = c - a * g / (2 * K + 1) + b * l / (2 * m + 1)

            k = math.sqrt(4 * coil1.r_avg * coil2.r_avg / (pow(coil1.r_avg + coil2.r_avg, 2) + pow(z, 2)))
            k_k = sp.special.ellipk(pow(k, 2))
            k_e = sp.special.ellipe(pow(k, 2))

            F_sum = F_sum + (inputs.mu0 * z * coil1.current * coil2.current * k / (4 * math.sqrt(coil1.r_avg * coil2.r_avg))
                             * ((2 - pow(k, 2)) / (1 - pow(k, 2)) * k_e - 2 * k_k))

    F_coeff = coil1.turns * coil2.turns / ((2 * K + 1) * (2 * m + 1))
    F = F_coeff * F_sum
    return F


def run_kinematics(projectile, F, timestep):
    acceleration = F / projectile.mass
    projectile.vel = projectile.vel + acceleration * timestep
    projectile.pos = projectile.pos + projectile.vel * timestep
