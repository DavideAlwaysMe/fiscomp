import matplotlib.pyplot as plt
import numpy as np

t, x, v, e, delta_E = np.loadtxt(
    'Runge_kutta_dt0.05.dat', unpack=True, usecols=(0, 1, 2, 3, 4))

# tempo e coordinata x
plt.plot(t, x, color='g')
plt.title('Algoritmo')
plt.xlabel('t')
plt.ylabel('x')
# plt.figure(figsize=(7,7))
plt.savefig('txx.png', dpi=199)
plt.clf()

# tempo e coordinata v
plt.plot(t, v, color='g')
plt.title('Algoritmo')
plt.xlabel('t')
plt.ylabel('v')
plt.savefig('txv.png', dpi=199)
plt.clf()

# x e v
plt.plot(x, v, color='g')
plt.title('Algoritmo')
plt.xlabel('x')
plt.ylabel('v')
plt.savefig('xxv.png', dpi=199)
plt.clf()
