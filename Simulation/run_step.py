import methods

on_drivers = []
driver_currents = []
def run_turn(timestep, stage_list, projectile):
    armature = projectile.coil
    armature.pos = projectile.pos
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
        methods.solve_energy(driver, timestep)

        force += methods.magnetic_force(stage, projectile)

    methods.increment_time(timestep, armature)
    methods.solve_thermals(armature, timestep)
    methods.solve_resistance(armature)
    methods.solve_kinematics(projectile, force, timestep)


    for stage in stage_list:
        driver = stage.coil
        driver.currents_list.append(driver.current)
        driver.energy_list.append(driver.energy)

    armature.energy_list.append(armature.energy)



    return projectile, armature
