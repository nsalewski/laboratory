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

params, covariance = curve_fit(f, w, phi)  # bounds = ([0], [np.inf]))
errors = np.sqrt(np.diag(covariance))
print('c =', params[0], '±', errors[0])

print(params)
ascii.write([w, a, b, phi], 'Messdaten/c.tex', format="latex")
ascii.write([np.round(phi, 5), np.round(urc / unull, 5)],
            'Messdaten/pol.tex', format='latex')
# used temp, bc without temp there was really freaky and wrong behaviour
# in matplotlib
m = np.logspace(0.01, 4)
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
phi_ = np.linspace(0, np.pi / 2, 1000)
v = -np.tan(phi_)
phi_ = np.linspace(0, np.pi / 2, 1000)
v = -np.tan(phi_) / (3.67 * 10**(-3))
phi_ = np.linspace(0.000000000000001, np.pi / 2, 1000)
v = -np.tan(phi_)
plt.polar(phi, urc / unull, 'rx', label='Messwerte')
plt.polar(phi_, -np.sin(phi_) / v, 'b-', label='Theoriekurve')
#plt.polar(phi_, -np.sin(phi_) / (v * 3.67 * 10 ** (-3)), 'b-', label = 'Theoriekurve')
#plt.polar(phi_, -np.sin(phi_) / (v * 3.67 * 10 ** (-3)), 'b-', label = 'Theoriekurve')
#plt.polar(phi_, -np.sin(phi_) / (v * 3.67 * 10 ** (-3)), 'b-', label = 'Theoriekurve')
xT = plt.xticks()[0]
xL = ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',
      r'$\pi$', r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$']
plt.xticks(xT, xL)
#([<matplotlib.axis.XTick object at 0x107bac490>, <matplotlib.axis.XTick object at 0x109a31310>, <matplotlib.axis.XTick object at 0x109a313d0>, <matplotlib.axis.XTick object at 0x109a31050>, <matplotlib.axis.XTick object at 0x1097a8690>, <matplotlib.axis.XTick object at 0x1097a8cd0>, <matplotlib.axis.XTick object at 0x1097a8150>, <matplotlib.axis.XTick object at 0x107bb8fd0>], <a list of 8 Text xticklabel objects>)
plt.tight_layout()
plt.polar(phi_, -np.sin(phi_) / v, 'b-', label='Theoriekurve')
plt.savefig('polaar.pdf')

print(v)
