import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii


w, urc, ug, a = np.genfromtxt("Messdaten/b_c.txt", unpack=True)
unull = 6.04
a = a / 1000
b = 1 / w
phi = 2 * np.pi * a / b


def f(w, c):
    return np.arctan(-w * c)

params, covariance = curve_fit(f, w, phi) # bounds = ([0], [np.inf]))
errors = np.sqrt(np.diag(covariance))
print('c =', params[0], '±', errors[0])
print(params)
ascii.write([w, a, b, phi], 'Messdaten/c.tex', format="latex")
ascii.write([phi, urc/unull], 'Messdaten/pol.tex', format = 'latex')
# used temp, bc without temp there was really freaky and wrong behaviour
# in matplotlib
m = np.logspace(0.01, 4,)
temp = (f(m, *params))
plt.plot(w, phi, 'rx', label="Messwerte")
plt.plot(m, temp, 'b-', label='Ausgleichskurve')
plt.xlim(4.24, 10000)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\phi(\omega)$/ $\si{\radian}$")
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/phase.pdf")

plt.clf()
phi_ = np.linspace(0, np.pi/2, 1000)
v = -np.tan(phi_) / (3.67 * 10**(-3))
plt.polar(phi, urc/unull, 'rx', label = 'Messwerte')
plt.polar(phi_, -np.sin(phi_) / (v * 3.67 * 10 ** (-3)), 'b-', label = 'Theoriekurve')
plt.savefig('polaar.pdf')
