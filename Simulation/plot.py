import matplotlib.pyplot as plt
import numpy as np
import main
import run_step

runtime = inputs.runtime
timestep = inputs.timestep

# figp, p = plt.subplots()
# p.plot(np.arange(0, runtime, timestep), positions)
# p.set(xlabel='Time (s)', ylabel='Position (m)', title='Projectile Position')
# figp.savefig("pos.png")

figv, v = plt.subplots()
v.plot(np.arange(0, runtime, timestep), velocities)
v.set(xlabel='Time (s)', ylabel='Velocity (m/s)', title='Projectile Velocity')
figv.savefig("vel.png")

figa, p = plt.subplots()
p.plot(np.arange(0, runtime, timestep), accelerations)
p.set(xlabel='Time (s)', ylabel='Acceleration (m/s^2)', title='Projectile Acceleration')
figa.savefig("accel.png")

figi, i = plt.subplots()
i.plot(np.arange(0, runtime, timestep), driver_current)
i.plot(np.arange(0, runtime, timestep), armature_current)
i.set(xlabel='Time (s)', ylabel='Current (A)', title='Driver and Armature Currents')
i.legend(['Driver', 'Armature'])
figi.savefig("current.png")

figth, th = plt.subplots()
th.plot(np.arange(0, runtime, timestep), driver_temp)
th.plot(np.arange(0, runtime, timestep), armature_temp)
th.set(xlabel='Time (s)', ylabel='Temperature (K)', title='Driver and Armature Temperatures')
th.legend(['Driver', 'Armature'])
figth.savefig("temperatures.png")

fige, e = plt.subplots()
e.plot(np.arange(0, runtime, timestep), driver_energy)
e.plot(np.arange(0, runtime, timestep), armature_energy)
e.set(xlabel='Time (s)', ylabel='Energy (J)', title='Energy Balance')
e.legend(['Coil Electrical Energy', 'Projectile Kinetic Energy'])
fige.savefig("energy.png")

figeff, eff = plt.subplots()
eff.plot(np.arange(0, runtime, timestep), efficiency)
eff.set(xlabel='Time (s)', ylabel='Efficiency', title='Energy Efficiency')
figeff.savefig("efficiency.png")

figm, m = plt.subplots()
m.plot(np.arange(0, runtime, timestep), m_inductances)
m.set(xlabel='Time (s)', ylabel='Mutual Inductance (H)', title='Mutual Inductance')
figm.savefig("m_inductance.png")
