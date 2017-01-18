import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
malu, dalu, dxalu = np.genfromtxt("Messdaten/a.txt", unpack=True)
D_x = dxalu - dalu
xalu = malu + 3
L = 0.54
m = 767.4 / 1000
g = 9.81  # scipy einfügen


def f(x, E):
    return ((L * x**2 - ((x**3) / 3)) * (F / (2 * I * E)))

Lstabalu = 60
dstabalu = 1
vstabalu = dstabalu**2 * Lstabalu
mstabalu = 167.1
pdichtealu = mstabalu / vstabalu
print("Dichte Stab rechteckig", pdichtealu)

Lstabmessing = 55
dstabmessing = 1
rstabmessing = dstabmessing / 2
mstabmessing = 360.5
vstabmessing = np.pi * rstabmessing**2 * Lstabmessing
pdichtemessing = mstabmessing / vstabmessing
print("Dichte Stab rund=", pdichtemessing)
# plt.plot()
