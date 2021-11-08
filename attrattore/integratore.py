import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

t, x, y, z, vx, vy, vz = np.loadtxt(
    'integratore.dat', usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, label='Traiettoria')
plt.savefig('traiettoria.png')

plt.clf()
raggio = np.sqrt(x**2+y**2+z**2)
plt.plot(t[:40000], raggio[:40000])
plt.savefig('raggio.png')

plt.clf()
plt.plot(t[40000:50000], x[40000:50000])
plt.grid()
plt.savefig('x.png')

plt.clf()
plt.plot(y, vy)
plt.show()
