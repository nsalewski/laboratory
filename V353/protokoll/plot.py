import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
#import sympy
#from uncertainties import ufloat
#import uncertainties.unumpy as unp
#from sympy import Symbol, latex
#from sympy import *
#from pylab import *

uc, t = np.genfromtxt("Messdaten/a.txt", unpack=True)
unull = 6.04
t = t / 1000
ascii.write([uc, np.round(np.log(uc / unull), 4), t], 'Messdaten/a.tex', format="latex",
            names=[r'$U_\\text{C}$/$\\si{\\volt}$', '$\\ln{(\\frac{U_\\text{C}}{U_\\text{0}})}$', "$t$ /$10^{-6}\\si{\\second}$"])


def f(x, a, b):
    return (-1 / a) * x + b
params, covariance = curve_fit(f, t, np.log(uc / unull))
errors = np.sqrt(np.diag(covariance))
print('a = ', params[0], '±', errors[0])
print('b=',params[1], '±', errors[1])
print(params)
m = np.linspace(-0.1, 2)
plt.plot(t, np.log(uc / unull), 'rx', label="Messwerte")
plt.plot(m, f(m, *params), 'b-', label='Ausgleichsgerade')
plt.xlabel(r'$t$ /$10^{-3}\si{\second}$')
plt.ylabel(r'$\ln\left(\frac{U_\text{C}}{U_\text{0}}\right)$')
plt.legend(loc='best')
plt.tight_layout()

plt.xlim([0, 2])
plt.savefig("build/test.pdf")
