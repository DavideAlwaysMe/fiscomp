import matplotlib.pyplot as plt
import numpy as np

# GRAFICO EULERO DT E DELTA_E CON FIT LINEARE

# PER OGNI ALGORITMO UN GRAFICO CON X E V
t, x, v, e, delta_E = np.loadtxt(
    'eulero_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Eulero')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('eulero.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'eulero-cromer_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Eulero-Cromer')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('eulero-cromer.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'mezzo_passo_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Mezzo Passo')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('mezzo_passo.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'punto_centrale_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Punto Centrale')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('punto_centrale.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'verlet_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Verlet')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('verlet.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'verlet_autosufficiente_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Verlet Autosufficiente')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('verlet_autosufficiente.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'predizione_correzione_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Predizione Correzione')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('predizione_correzione.png', dpi=199)
plt.clf()

t, x, v, e, delta_E = np.loadtxt(
    'runge_kutta_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 3, 4))
plt.title('Runge Kutta')
plt.xlabel('t')
plt.ylabel('x(t) e v(t)')
# tempo e coordinata x
plt.plot(t, x, color='b', label='x')
# tempo e coordinata v
plt.plot(t, v, color='c', label='v')
plt.legend()
plt.savefig('runge_kutta.png', dpi=199)
plt.clf()

# GRAFICO X E V EULERO ED EULERO CROMER PER VERIFICARE LA DIVERGENZA DI EULERO
# x e v
t1, x1, v1, delta_E1 = np.loadtxt(
    'eulero_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 4))
t2, x2, v2, delta_E2 = np.loadtxt(
    'eulero-cromer_dt0.01.dat', unpack=True, usecols=(0, 1, 2, 4))
plt.figure(figsize=(7, 7))
plt.plot(x1, v1, color='orange', label='Eulero')
plt.plot(x2, v2, color='red', label='Eulero-Cromer')
plt.title('Divergenza di Eulero ed Eulero Cromer')
plt.xlabel('x')
plt.ylabel('v')
plt.legend()
plt.savefig('divergenza_eulero.png', dpi=199)
plt.clf()

# CONFRONTO EULERO ED EULERO CROMER CON LE X (CAMBIARE T O DT PER NOTARE ERRORE)
plt.figure(figsize=(7, 7))
plt.plot(t1, x1, color='orange', label='Eulero')
plt.plot(t2, x2, color='red', label='Eulero-Cromer')
plt.title('Confronto tra Eulero ed Eulero Cromer')
plt.xlabel('x')
plt.ylabel('v')
plt.legend()
plt.savefig('eulero_vs_eulero-cromer.png', dpi=199)
plt.clf()

# GRAFICO EULERO CROMER TEMPO E DELTA
plt.figure(figsize=(7, 7))
plt.plot(t1, delta_E1, color='orange', label='Eulero')
plt.plot(t2, delta_E2, color='red', label='Eulero-Cromer')
plt.title('Confronto tra Eulero ed Eulero Cromer')
plt.xlabel('t')
plt.ylabel('$\Delta$E(t)/E(0)')
plt.legend()
plt.savefig('eulero_eulero-cromer_deltaE.png', dpi=199)
plt.clf()

# GRAFICO EULERO CROMER TEMPO E DELTA ZOOM PER VERIFICARE CHE OSCILLA
plt.figure(figsize=(7, 7))
plt.plot(t1, delta_E1, color='orange', label='Eulero')
plt.plot(t2, delta_E2, color='red', label='Eulero-Cromer')
plt.title('Confronto tra Eulero ed Eulero Cromer')
plt.xlabel('t')
plt.ylabel('$\Delta$E(t)/E(0)')
plt.xlim([0, 4])
plt.ylim([-0.02, 0.02])
plt.legend()
plt.savefig('eulero_eulero-cromer_deltaE_zoom.png', dpi=199)
plt.clf()
