import methods

on_drivers = []
driver_currents = []
def run_turn(timestep, stage_list, projectile):
    armature = projectile.coil
    armature.pos = projectile.pos
    armature.vel = projectile.vel
    switch = min(projectile.stage, len(stage_list) - 1)
    force = 0

    if(armature.pos >= stage_list[switch].trigger_pos and projectile.stage < len(stage_list)):
        on_drivers.append(stage_list[switch])
        projectile.stage += 1
        print("Switched stage", switch + 1)

    for stage in on_drivers:
        driver = stage.coil

        methods.increment_time(timestep, driver)
        methods.solve_currents(stage, projectile, driver.time_on, armature.time_on)
        methods.solve_thermals(driver, timestep)
        methods.solve_energy(driver, timestep)

        force += methods.magnetic_force(stage, projectile)

    methods.increment_time(timestep, armature)
    methods.solve_thermals(armature, timestep)
    methods.solve_resistance(armature)

    armature.currents_list.append(armature.current)
    armature.voltage_list.append(armature.voltage)
    armature.energy_list.append(armature.energy)
    armature.temperature_list.append(armature.temp)

    methods.solve_kinematics(projectile, force, timestep)

    for stage in stage_list:
        driver = stage.coil
        driver.currents_list.append(driver.current)
        driver.voltage_list.append(driver.voltage)
        driver.energy_list.append(driver.energy)
        driver.temperature_list.append(driver.temp)

    return projectile, armature
