import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
g = ufloat(9.811899, 0.000041)

x, d0, d = np.genfromtxt("Messdaten/c.txt", unpack=True)
D = d - d0
x1 = x[0:26]
x2 = x[28:52]
D1 = D[0:26]
D2 = D[28:52]

ascii.write([x, d0, d, D], 'Messdaten/beidseitig.tex', format="latex",
            names=["messpunkt x ", "D0", "Dlast", "D"])

m_alu, d_alu, dx_alu = np.genfromtxt("Messdaten/a.txt", unpack=True)
L_stab_alueingespannt = 0.56
mlast_alu = 4720.1
L_stab_alu = 60
d_stab_alu = 0.01
v_stab_alu = d_stab_alu**2 * L_stab_alu
m_stab_alu = 167.1
pdichte_alu = m_stab_alu / v_stab_alu
print("Dichte Stab rechteckig", pdichte_alu)
dichte_lit_alu = 2.712  # in g/m^3
print("Dichte Alu Literatur", dichte_lit_alu)


x = x / 100
x1 = x1 / 100
x2 = x2 / 100
x_alu_fit = ((3 * L_stab_alueingespannt**2) * x1 - 4 * x1**3)


def Y1(x, a):
    return a * x

params, covariance = curve_fit(Y1, x_alu_fit, D1)
errors = np.sqrt(np.diag(covariance))
print("params", *params, "und +/-", errors[0])
plt.plot(x_alu_fit, D1, 'rx', label="Messwerte")
plt.plot(x_alu_fit, Y1(x_alu_fit, *params), 'b-', label="Regressionsgrade")
plt.xlabel(
    r"$3L^2 x - 4x^3$")
plt.ylabel(
    r"$D(x)$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/c.pdf')
a_alu = ufloat(params[0], errors[0])
F_alu = mlast_alu * g
I_alu = d_stab_alu**4 / 12

E_alu = F_alu / (48 * a_alu * I_alu)
print("E alu=", E_alu)
##########################################################################


def Y2(x, A):
    return A * x

x_alufit = 4 * x2**3 - 12 * L_stab_alueingespannt * x2**2 + \
    9 * L_stab_alueingespannt**2 * x2 - L_stab_alueingespannt**3
plt.clf()
params, covariance = curve_fit(Y2, x_alufit, D2)
errors = np.sqrt(np.diag(covariance))
print("params", *params, "fehler", *errors)
plt.plot(x_alufit, D2, 'rx', label="Messwerte")
plt.plot(x_alufit, Y2(x_alufit, *params), 'b-', label="Regressionsgrade")
plt.xlabel(
    r"$4x^3 -12Lx^2 + 9L^2x - L^3")
plt.ylabel(
    r"$D(x)$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/c2.pdf')
a_alu = ufloat(params[0], errors[0])
F_alu = mlast_alu * g
I_alu = d_stab_alu**4 / 12

E_alu = F_alu / (48 * a_alu * I_alu)
print("E alu=", E_alu)
