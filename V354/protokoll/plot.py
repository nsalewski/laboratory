import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import sympy
from uncertainties import ufloat
import uncertainties.unumpy as unp
from sympy import Symbol, latex
from sympy import *
from pylab import *
uc, t = np.genfromtxt("Messdaten/a_and_b.txt", unpack=True)

ascii.write([uc, t], 'Messdaten/tab_a.tex', format="latex",
            names=[r'$U_\\text{C}$/$\\si{\\volt}$', "$t$ /$10^{-6}\\si{\\second}$"])

def f(x, a, b):
    return a*np.exp(-2*np.pi*b*t*(1/1000000))
params, covariance = curve_fit(f, t, uc)
errors = np.sqrt(np.diag(covariance))
print(errors)
print('A_0 = ', params[0], '±', errors[0])
print('mü = ', params[1], '±', errors[1])
plt.plot(t, uc, 'rx', label="Messwerte")
plt.plot(t, f(t, *params), 'b-', label='Ausgleichsgerade') #didnt work with linspace for me, dunno why.
plt.xlabel(r'$t$ /$10^{-6}\si{\second}$')
plt.ylabel(r'$U_\text{C}$/$\si{\volt}$')
plt.legend(loc='best')
plt.tight_layout()

plt.savefig("build/test.pdf")
