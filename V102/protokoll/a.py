import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants

# given values

m_k = ufloat(512.2, 0.2)
m_k = m_k / 1000
D_k = ufloat(50.760, 0.004)
D_k = D_k / 1000
R_k = D_k / 2
Theta_halterung = 22.5
Theta_halterung = Theta_halterung / 10**7
L_draht = 0.585
E = ufloat(210.0, 0.5) * 10**9

# calculation :moment of inertia
print("Teil a:")

Theta_kugel = 2 / 5 * m_k * R_k**2
print("Theta Kugel=", Theta_kugel)
Theta_ges = Theta_kugel + Theta_halterung
print("Theta_ges", Theta_ges)

# caltulation: R_draht

n, Ddraht = np.genfromtxt("Messdaten/daten_des_aufbau.txt", unpack=True)
Rdraht = Ddraht / 2
R_draht = ufloat(np.mean(Rdraht), np.std(
    Rdraht, ddof=1) / np.sqrt(len(Rdraht)))
R_draht = R_draht / 1000
ascii.write([n, Ddraht, Rdraht], 'Messdaten/draht.tex', format="latex",
            names=["messpunkt", "ddraht", "rdraht"])
print("R_draht = ", R_draht)

# calculation Period
m, t = np.genfromtxt("Messdaten/a.txt", unpack=True)
T = ufloat(np.mean(t), np.std(t, ddof=1) / np.sqrt(len(t)))
print("T=", T)
ascii.write([m, t], 'Messdaten/schwingungsdauer.tex', format="latex",
            names=["Messpunkt", "Schwingungsdauer"])
# calculation :shear modulus
G = ((2 * L_draht) / (np.pi * R_draht**4)) * (2 * np.pi / T)**2 * Theta_ges
print("Schubmodul G=", G)
# calculation elastic constants
muh = ((E / (2 * G)) - 1)
Q = (E) / (3 * (1 - 2 * muh))
print("E=", E)
print("müh=", muh)
print("Q=", Q)
print("***********************Ende von a******************************")

# calculation: part b

N = 390
R_Helmholtz = 78 / 1000

# calculation helmholtz-coils magnetic field
amperage = np.linspace(0.1, 1, 10)
mu_null = scipy.constants.physical_constants["mag. constant"]
mu_0 = mu_null[0]  # muh contains (value, "unit",error)
B_H = 8 / np.sqrt(125) * amperage * N * mu_0 / R_Helmholtz

###########
S, Ts = np.genfromtxt("Messdaten/b.txt", unpack="True")
mean = list()
t1 = list()  # these lists are just for proper formatting and writing tables
t2 = list()
t3 = list()
t4 = list()
t5 = list()
for i in range(10):
    mean.append(ufloat(np.mean(Ts[(i * 5): ((i * 5) + 5)]), np.std(
        Ts[(i * 5): ((i * 5) + 5)], ddof=1) / np.sqrt(len(Ts[(i * 5): ((i * 5) + 5)]))))
    t1.append(Ts[(i * 5)])
    t2.append(Ts[((i * 5) + 1)])
    t3.append(Ts[((i * 5) + 2)])
    t4.append(Ts[((i * 5) + 3)])
    t5.append(Ts[((i * 5) + 4)])
# array with averages of T for each amperage
mean_T_Per_Amperage = np.asarray(mean)


# right side of equation (formula (20) => m*B+D=Phi) ;contains errors
# (Theta_ges and mean_T_Per_Amperage are uarrays)
# for plotting we will only need nominal values of Phi

# calculate regression and plotting for m*B+D=Phi

print(mean_T_Per_Amperage)
mean_T_Per_Amperage_nv = list()
for i in range(10):
    mean_T_Per_Amperage_nv.append(mean_T_Per_Amperage[i].nominal_value)
mean_T_Per_Amperage_nominalvalue = np.asarray(mean_T_Per_Amperage_nv)


def ausgleichsgrade(B, a, b):
    return a * B + b

params, covariance = curve_fit(
    ausgleichsgrade, (1 / mean_T_Per_Amperage_nominalvalue**2), B_H)
errors = np.sqrt(np.diag(covariance))
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
m = ((4 * np.pi**2) / a) * Theta_ges
D = -m * b
print('m =', m)
print('D = ', D)
ausgleich_T = np.linspace(5, 60, 100)
plt.plot(1 / mean_T_Per_Amperage_nominalvalue **
         2, B_H * 1000, 'ro', label="Messwerte")
plt.plot(1 / ausgleich_T**2, ausgleichsgrade(
    1 / ausgleich_T**2, *params) * 1000, 'b-', label="Regressionsgrade")
plt.xlabel(
    r"$\frac{1}{T_{\mathrm{m}}^2}$/$\si{\per\square\second}$")
plt.ylabel(
    r"$B_{\mathrm{Helmholtz}} \cdot 10^{-3}$/ $\si{\tesla}$")
plt.ylim(0, 5)
plt.xlim(0, 0.035)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b.pdf')

# generating tables
ascii.write([amperage, t1, t2, t3, t4, t5, mean_T_Per_Amperage, B_H], 'Messdaten/b.tex', format="latex",
            names=["Stromstärke", "T1", "T2", "T3", "T4", "T5", "T mittel \pm standardabweichung", "BFeld"])
