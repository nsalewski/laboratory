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

i,u=np.genfromtxt("Messdaten/b.txt", unpack=True)
ascii.write([u, i], 'Messdaten/b.dat', format="latex", names=['$U_\text{k}$/$\si{\volt}$', '$I$/$10^{-3}\si{\ampere}$'])
def f(x, a, b):
    return a * x + b

params, covariance = curve_fit(f,i, u)
errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
x=np.linspace(20,100)

plt.plot(i, u, 'rx', label="Messwerte")
plt.plot(x, f(x, *params), 'b-', label='Linearer Fit')
plt.xlabel("Stromstärke $I$/$\si{\\milli\\ampere}$")
plt.ylabel("Klemmenspannung $U_\\text{k}$/$\si{\\milli\\volt}$")
plt.legend(loc="best")
plt.savefig("Bilder/a.pdf")
