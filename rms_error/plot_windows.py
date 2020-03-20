import numpy as np
import matplotlib.pyplot as plt
import time
import dask.dataframe as dd

'''
Identifier:
    c: cosine coefficients
    s: sine coefficients
    fs: sampling frequency
    sec_whole_day: seconds in a whole day
    n: sampling numbers
    position_d: gracefo d satellite positions in a whole day 
    velocity_d: gracefo d satellite velocities in a whole day
    gracefo_rpt: the data from integrator programme
'''

time_start = time.time()

# sampling interval = 5s
# for now, data of a whole day is considered
fs = 10

gracefo_rpt = np.array(
    dd.read_csv('..//output//report.txt', sep='\s+', usecols=[0, 1, 2, 3, 4, 5, 6], dtype=np.float64)
)

n, temp = np.shape(gracefo_rpt)
sec_whole_day = int(n * fs)

position_velocity_d = np.loadtxt('gracefo-d-2019-0301.txt').astype(np.float64)

grav = np.array(
    dd.read_csv('..//output//potential.txt', dtype=np.float64)
)

# grace fo true position vector and velocity vector
index_span = np.linspace(0, sec_whole_day - fs, n, dtype=np.int)
position_d = position_velocity_d[index_span, 1: 4]
velocity_d = position_velocity_d[index_span, 4: 8]

# whether the energy is balanced or not
delta_grav = grav[1: n] - grav[0]
delta_kine = (
                     - np.square(np.linalg.norm(gracefo_rpt[1: n, 4: 8], axis=1))
                     + np.square(np.linalg.norm(gracefo_rpt[0,  4: 8]))
             ) / 2.
delta_energy = delta_kine + delta_grav[:, 0]
np.savetxt('..//output//kine_energy.txt', np.square(np.linalg.norm(gracefo_rpt[:, 4: 8], axis=1)))
np.savetxt('delta_energy_5s.txt', delta_energy, newline='\n')

plt.figure('fig_energy')
plt.plot(index_span[1: n], delta_energy, linewidth=1)
plt.xlabel('Time (s)')
plt.ylabel('Energy difference')
plt.title('Energy difference between the $t_n$ epoch and $t_0$ epoch', fontsize=10)
plt.show()

# calculate the rms error of velocity and position
# in the reference, the chosen time span is about 90min, aka, 5400s.
points = n
orbits = np.floor(sec_whole_day / 5400.)

rms_position = np.linalg.norm(
    (np.linalg.norm(position_d[0: points], axis=1) - np.linalg.norm(gracefo_rpt[0: points, 1: 4], axis=1))**2 / points
)
print(np.linalg.norm(position_d[0: points], axis=1) - np.linalg.norm(gracefo_rpt[0: points, 1: 4], axis=1))
error_ratio_position = rms_position / orbits / np.max(np.linalg.norm(position_d))

rms_velocity = np.linalg.norm(
    (np.linalg.norm(velocity_d[0: points], axis=1) - np.linalg.norm(gracefo_rpt[0: points, 4: 8], axis=1))**2 / points
)
error_ratio_velocity = rms_velocity / orbits / np.max(np.linalg.norm(velocity_d))

print('The root mean square of the position error ratio is ', error_ratio_position)
print('The root mean square of the velocity error ratio is ', error_ratio_velocity)

print('This programme costs ', time.time() - time_start, '.')
