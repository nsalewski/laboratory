import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants
m_k = ufloat(512.2, 0.2)
m_k = m_k / 1000
D_k = ufloat(50.760, 0.004)
D_k = D_k / 1000
R_k = D_k / 2
Theta_halterung = 22.5
Theta_halterung = Theta_halterung / 10**7
L_draht = 0.585
Theta_kugel = 2 / 5 * m_k * R_k**2
E = ufloat(210.0, 0.5) * 10**9
print("Teil a:")
print("Theta Kugel=", Theta_kugel)
Theta_ges = Theta_kugel + Theta_halterung
n, Ddraht = np.genfromtxt("Messdaten/daten_des_aufbau.txt", unpack=True)
Rdraht = Ddraht / 2
R_draht = np.mean(Rdraht)
R_draht = R_draht / 1000
m, t = np.genfromtxt("Messdaten/a.txt", unpack=True)
T = np.mean(t)

G = ((2 * L_draht) / (np.pi * R_draht**4)) * (2 * np.pi / T)**2 * Theta_ges
print("Schubmodul G=", G)
muh = ((E / (2 * G)) - 1)
Q = (E) / (3 * (1 - 2 * muh))
print("E=", E)
print("müh=", muh)
print("Q=", Q)
print("***********************Ende von a******************************")
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
mu_0 = muh[0]  # muh contains (value, "unit",error)
B_H = 8 / np.sqrt(125) * amperage * N * mu_0 / R_Helmholtz
Delta = ((2 * np.pi) / mws)**2 * Theta_ges.nominal_value
print("B=", B_H)
print("Theta_ges ohne Fehler", Theta_ges.nominal_value)


def ausgleichsgrade(B, a, b):
    return a * B + b


params, covariance = curve_fit(ausgleichsgrade, B_H, Delta)

errors = np.sqrt(np.diag(covariance))


m = ufloat(params[0], errors[0])
D = ufloat(params[1], errors[1])
print('a = m =', params[0], '+-', errors[0])
print('b = D = ', params[1], '+-', errors[1])

ascii.write([amperage, t1, t2, t3, t4, t5, mws], 'Messdaten/b.tex', format="latex",
            names=["Stromstärke", "T1", "T2", "T3", "T4", "T5", "T mittel"])

ausgleich_B = np.linspace(0, 50 / 10000, 100)
plt.plot(B_H * 1000, Delta * 10000, 'ro', label="Messwerte")
plt.plot(ausgleich_B * 1000, ausgleichsgrade(
    ausgleich_B, *params) * 10000, 'b-', label="Theoriekurve")
plt.ylabel(r"$Vielkram*10^{-5}$")
plt.xlabel(r"$B*10^{-3}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b.pdf')
