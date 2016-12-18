import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii

w, a, uc, u = np.genfromtxt("Messdaten/c_und_d.txt", unpack=True)
b = 1 / (w * 1000)
phi = 2 * np.pi * a / (10**6 * b)
R = 559.5 #müsste doch 559.5 sein, wegen generatorinnenwiderstand. habs mal geändert. Vorher war es 509.5
L = 10.11 / 1000
C = 2.098 / 10**9

def f(w):
    return np.arctan(-(w * R * C) / (1 - L * C * w**2))

ascii.write([w, a, b, phi], 'Messdaten/d.tex', format="latex")

#m = np.logspace(0.01, 4)
#temp = (f(m, *params))

plt.figure(0)
plt.plot(w, phi, 'rx', label="Messwerte")
#plt.plot(m, temp, 'b-', label='Ausgleichskurve')
plt.xlim(1, 100)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\phi(\omega)$/ $\si{\radian}$")
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/taskd.pdf")

plt.figure(1)
plt.plot(w, phi, 'rx', label="Messwerte")
plt.axhline(y = np.pi / 4)
plt.axhline(y = np.pi / 2)
plt.axhline(y = 3 * np.pi / 4)
plt.xlim(10, 40)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\phi(\omega)$/ $\si{\radian}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/taskdlinear.pdf")

