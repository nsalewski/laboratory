import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants

N = 390
R_Helmholtz = 78 / 1000


###########
S, Ts = np.genfromtxt("Messdaten/b.txt", unpack="True")
mean = list()
t1 = list()
t2 = list()
t3 = list()
t4 = list()
t5 = list()
for i in range(10):
    mean.append(np.mean(Ts[(i * 5): ((i * 5) + 5)]))
    t1.append(Ts[(i * 5)])
    t2.append(Ts[((i * 5) + 1)])
    t3.append(Ts[((i * 5) + 2)])
    t4.append(Ts[((i * 5) + 3)])
    t5.append(Ts[((i * 5) + 4)])
mws = np.asarray(mean)  # array with averages of T for each amperage
amperage = np.linspace(0.1, 1, 10)
muh = scipy.constants.physical_constants["mag. constant"]
mu_0 = ufloat(muh[0], muh[2])  # muh contains (value, "unit",error)
B_H = 8 / np.sqrt(125) * amperage * N * mu_0 / R_Helmholtz


def ausgleichsgrade(B, a, b):
    return a * B + b


params, covariance = curve_fit(B_H, Delta, ausgleichsgrade)

errors = np.sqrt(np.diag(covariance))

m = ufloat(params[0], errors[0])
D = ufloat(params[1], errors[1])
print('a = m =', params[0], '+-', errors[0])
print('b = D = ', params[1], '+-', errors[1])

ascii.write([amperage, t1, t2, t3, t4, t5, mws], 'Messdaten/b.tex', format="latex",
            names=["Stromstärke", "T1", "T2", "T3", "T4", "T5", "T mittel"])

ausgleich_B = np.linspace(0, 93, 100)
plt.plot(B_H, Delta, 'ro', label="Messwerte")
plt.plot(ausgleich_B, ausgleichsgrade(
    ausgleich_B, *params), 'b-', label="Theoriekurve")
plt.ylabel(r"$Vielkram$")
plt.xlabel(r"$B$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b.pdf')
