import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

twobody_rpt = np.loadtxt('report.txt')
coor = np.loadtxt('coordinates_earth_fixed_system.txt', delimiter='\t').astype(np.float64)

[raw, col] = np.shape(twobody_rpt)

r = np.zeros(raw)
v = np.zeros(raw)
error = np.zeros(raw)
mu = 398600.4415

for index in range(raw):
    r[index] = np.linalg.norm([twobody_rpt[index, 1], twobody_rpt[index, 2], twobody_rpt[index, 3]])
    v[index] = np.linalg.norm([twobody_rpt[index, 4], twobody_rpt[index, 5], twobody_rpt[index, 6]])
    error[index] = (v[index]**2 - v[0]**2 + 2 * mu / r[0] - 2 * mu / r[index])

plt.figure()
plt.plot(twobody_rpt[:, 0], error)
plt.xlabel('time/s', fontsize=20)
plt.ylabel('error', fontsize=20)
plt.title('Spectral error (interval=100s) ($v_t^2-v_0^2 + 2GM(1/r_0-1/r_t))$', fontsize=20)
# plt.savefig('\\image\\error.png')
plt.show()

plt.figure()
ax_3 = plt.axes(projection='3d')
ax_3.plot3D(twobody_rpt[:, 1], twobody_rpt[:, 2], twobody_rpt[:, 3])
plt.show()
