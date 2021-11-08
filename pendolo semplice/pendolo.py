import matplotlib.pyplot as plt
import numpy as np

# GRAFICO EULERO DT E DELTA_E CON FIT LINEARE

# PER OGNI ALGORITMO UN GRAFICO CON X E V
t, x, v, e, delta_E = np.loadtxt(
    'pendolosemplice.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Eulero')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('ciao.png', dpi=199)
plt.clf()
