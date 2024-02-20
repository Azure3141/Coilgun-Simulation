import methods

on_drivers = []
def run_turn(timestep, stage_list, projectile):

    switch = min(projectile.stage, len(stage_list) - 1)
    force = 0

    if(projectile.pos >= stage_list[switch].trigger_pos and switch <= len(stage_list)):
        on_drivers.append(stage_list[switch])
        projectile.stage += 1
        print("Switched stage", switch + 1)

    for stage in on_drivers:
        driver = stage.coil
        armature = projectile.coil

        methods.increment_time(timestep, driver)
        methods.increment_time(timestep, armature)
        methods.solve_currents(stage, projectile, driver.time_on, armature.time_on)
        methods.solve_thermals(driver, timestep)
        methods.solve_thermals(armature, timestep)
        methods.solve_resistance(armature)

        force += methods.magnetic_force(stage, projectile)

    methods.solve_kinematics(projectile, force, timestep)
    #print(stage_list[0].coil.current, stage_list[1].coil.current)
    return projectile
