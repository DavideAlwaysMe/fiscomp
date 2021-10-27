import matplotlib.pyplot as plt
import numpy as np

# PER OGNI ALGORITMO UN GRAFICO CON X E V
# CONFRONTO EULERO ED EULERO CROMER CON LE X (CAMBIARE T O DT PER NOTARE ERRORE)
# GRAFICO X E V EULERO ED EULERO CROMER PER VERIFICARE LA DIVERGENZA DI EULERO
# GRAFICO EULERO TEMPO ED DELTA E PER VERIFICARE LA LINEARITÀ
# GRAFICO EULERO CROMER TEMPO E DELTA E PER VERIFICARE CHE OSCILLA
# GRAFICO EULERO DT E DELTA_E CON FIT LINEARE

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
