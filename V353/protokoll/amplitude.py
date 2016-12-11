import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii


w, urc, ug, phi = np.genfromtxt("Messdaten/b_c.txt", unpack=True)
unull = 7.28


def f(x, a):
    return (1 / (np.sqrt(1 + (x**2) * (a**2))))
params, covariance = curve_fit(f, w, urc / unull)
errors = np.sqrt(np.diag(covariance))
print('a = (RC) =', params[0], '±', errors[0])
ascii.write([urc, np.round(urc / unull, 2), w],
            'Messdaten/b.tex', format="latex")

# used temp, bc without temp there was really freaky and wrong behaviour
# in matplotlib
m = np.logspace(0.01, 4,)
temp = (f(m, *params))
plt.plot(w, urc / unull, 'rx', label="Messwerte")
plt.plot(m, temp, 'b-', label='Regression')
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\frac{A(\omega)}{U_\text{0}}$")
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/amplitude.pdf")
