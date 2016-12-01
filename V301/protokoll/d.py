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
i,u=np.genfromtxt("Messdaten/rechteckspannung.txt", unpack=True)
u=u*1000

def f(x, a, b):
    return a * x + b

params, covariance = curve_fit(f,i, u)
errors = np.sqrt(np.diag(covariance))


x=np.linspace(0,5)
plt.plot(i, u, 'rx', label="Messwerte")
plt.plot(x, f(x, *params), 'b-', label='Linearer Fit')
plt.xlabel("Stromstärke $I$/$\\si{\\milli\\ampere}$")
plt.ylabel("Klemmenspannung $U_\\text{k}$/$\si{\\milli\\volt}$")
plt.legend(loc="best")
plt.savefig("Bilder/d.pdf")
