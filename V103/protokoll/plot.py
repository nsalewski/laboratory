import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
m_alu, d_alu, dx_alu = np.genfromtxt("Messdaten/a.txt", unpack=True)
m_messing, d_messing, dx_messing = np.genfromtxt(
    "Messdaten/b.txt", unpack=True)
D_x_messing = dx_messing - d_messing
x_messing = m_messing + 3
mlast_messing = 520.9

D_x_alu = dx_alu - d_alu
x_alu = m_alu + 3
L_stab_alueingespannt = 0.54
mlast_alu = 767.4
# g = 9.81  # scipy einfügen

L_stab_alu = 60
d_stab_alu = 1
v_stab_alu = d_stab_alu**2 * L_stab_alu
m_stab_alu = 167.1
pdichte_alu = m_stab_alu / v_stab_alu
print("Dichte Stab rechteckig", pdichte_alu)
dichte_lit_alu = 2.712  # in g/m^3
print("Dichte Alu Literatur", dichte_lit_alu)

L_stab_messing = 55
d_stab_messing = 1
r_stab_messing = d_stab_messing / 2
m_stab_messing = 360.5
L_stab_messing_eingespannt = 0.49
v_stab_messing = np.pi * r_stab_messing**2 * L_stab_messing
pdichte_messing = m_stab_messing / v_stab_messing
print("Dichte stab rund=", pdichte_messing)
dichte_lit_messing = 8.400  # jeweils nach engineers toolbox
print("Dichte Messing Literatur", dichte_lit_messing)


ascii.write([x_alu, d_alu, dx_alu, D_x_alu], 'Messdaten/alu_einseitig.tex', format="latex",
            names=["messpunkt x ", "D_0", "D_belastet", "D diff"])

# Bisher: Plot von D(x). Soll aber plot von D(x_alu_fit) sein. Warum?!

x_alu = x_alu / 100
x_messing = x_messing / 100
x_alu_fit = (L_stab_alueingespannt * x_alu**2 - ((x_alu**3) / 3))
x_messing_fit = (L_stab_messing_eingespannt *
                 x_messing**2 - ((x_messing**3) / 3))


def D(x, a, b):
    return a * x + b

params, covariance = curve_fit(D, x_alu_fit, D_x_alu)
errors = np.sqrt(np.diag(covariance))
print("params", *params)
plt.plot(x_alu_fit, D_x_alu, 'rx', label="Messwerte")
plt.plot(x_alu_fit, D(x_alu_fit, *params), 'b-', label="Regressionsgrade")
plt.xlabel(
    r"$L\cdot x^2- \frac{x^3}{3}$")
plt.ylabel(
    r"$D(x)$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/a.pdf')
plt.clf()
params_m, covariance_m = curve_fit(D, x_messing_fit, D_x_messing)
errors = np.sqrt(np.diag(covariance_m))
print("params", *params_m)
plt.plot(x_messing_fit, D_x_messing, 'rx', label="Messwerte")
plt.plot(x_messing_fit, D(x_messing_fit, *params_m),
         'b-', label="Regressionsgrade")
plt.xlabel(
    r"$L\cdot x^2- \frac{x^3}{3}$")
plt.ylabel(
    r"$D(x)$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b.pdf')
