import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


N = 390
R_Helmholtz = 78 / 1000
#B_H = 8 / np.sqrt(125) * I * N * muhnull / R_Helmholtz


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


# def theorie(theta):
#    return np.sqrt(2 / (L * C) * (1 - np.cos(theta)))

ascii.write([amperage, t1, t2, t3, t4, t5, mws], 'Messdaten/b.tex', format="latex",
            names=["Stromstärke", "T1", "T2", "T3", "T4", "T5", "T mittel"])


#plt.plot(theta, w / 1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot) / 1000, 'b-', label="Theoriekurve")
#
# plt.ylabel(r"$Vielkram$")
# plt.xlabel(r"$$")
# plt.legend(loc='best')
# plt.tight_layout()
# plt.savefig('Bilder/b.pdf')
