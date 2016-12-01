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
i2,u2=np.genfromtxt("Messdaten/c.txt", unpack=True)
ir,ur=np.genfromtxt("Messdaten/rechteckspannung.txt", unpack=True)
i_s,us=np.genfromtxt("Messdaten/sinusspannung.txt", unpack=True)
u2=u2*1000
ur=ur*1000
us=us*1000
ascii.write([u, i], 'Messdaten/b.dat', format="latex", names=['$U_\text{k}$/$\si{\volt}$', '$I$/$10^{-3}\si{\ampere}$'])
ascii.write([u2, i2], 'Messdaten/c.tex', format="latex", names=['$U_\\text{k}$/$\\si{\\milli\\volt}$', '$I$/$\\si{\\milli\\ampere}$'])
ascii.write([ur, ir], 'Messdaten/r.tex', format="latex", names=['$U_\\text{k}$/$\\si{\\milli\\volt}$', '$I$/$\\si{\\milli\\ampere}$'])
ascii.write([us, i_s], 'Messdaten/s.tex', format="latex", names=['$U_\\text{k}$/$\\si{\\milli\\volt}$', '$I$/$\\si{\\milli\\ampere}$'])


def f(x, a, b):
    return a * x + b
def a(x,y):
    return(((np.mean(x*y)-np.mean(x)*np.mean(y))/(np.mean(x**2)-np.mean(x)**2)))
paramsb, covarianceb = curve_fit(f,i, u)
errorsb = np.sqrt(np.diag(covarianceb))
print('a Monozelle =', paramsb[0], '±', errorsb[0])
print('b Monozelle =', paramsb[1], '±', errorsb[1])
alpa=a(i_s,us)
print("test=", alpa)
paramsc, covariancec = curve_fit(f,i2, u2)
errorsc = np.sqrt(np.diag(covariancec))
print('a Monozelle 2=', paramsc[0], '±', errorsc[0])
print('b Monozelle 2=', paramsc[1], '±', errorsc[1])

paramsr, covariancer = curve_fit(f,ir, ur)
errorsr = np.sqrt(np.diag(covariancer))
print('a r=', paramsr[0], '±', errorsr[0])
print('b r=', paramsr[1], '±', errorsr[1])

paramss, covariances = curve_fit(f,i_s, us)
errorss = np.sqrt(np.diag(covariances))
print('a s=', paramss[0], '±', errorss[0])
print('b s=', paramss[1], '±', errorss[1])

x=np.linspace(20,100)
plt.plot(i, u, 'rx', label="Messwerte")
plt.plot(x, f(x, *paramsb), 'b-', label='Linearer Fit')
plt.xlabel("Stromstärke $I$/$\\si{\\milli\\ampere}$")
plt.ylabel("Klemmenspannung $U_\\text{k}$/$\si{\\milli\\volt}$")
plt.legend(loc="best")
plt.savefig("Bilder/a.pdf")
