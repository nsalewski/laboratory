import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
g = ufloat(9.811899, 0.000041)
x_linspace = np.linspace(0, 50) / 100


def D(x, a):
    return a * x


m_messing, d_messing, dx_messing = np.genfromtxt(
    "Messdaten/b.txt", unpack=True)
D_x_messing = dx_messing - d_messing
x_messing = m_messing + 3
mlast_messing = 520.9
L_stab_messing = 55
d_stab_messing = 0.01
r_stab_messing = d_stab_messing / 2
m_stab_messing = 360.5
L_stab_messing_eingespannt = 0.49
v_stab_messing = np.pi * r_stab_messing**2 * L_stab_messing
pdichte_messing = m_stab_messing / v_stab_messing
print("Dichte stab rund=", pdichte_messing)
dichte_lit_messing = 8.400  # jeweils nach engineers toolbox
print("Dichte Messing Literatur", dichte_lit_messing)
I_messing = np.pi * r_stab_messing**4 / 4
F_messing = mlast_messing * g
print("I messing=", I_messing)

x_messing = x_messing / 100
x_messing_fit = (L_stab_messing_eingespannt *
                 x_messing**2 - ((x_messing**3) / 3))
x_linspace_fit_messing = (L_stab_messing_eingespannt *
                          x_linspace**2 - ((x_linspace**3) / 3))
params_m, covariance_m = curve_fit(D, x_messing_fit, D_x_messing)
errors_m = np.sqrt(np.diag(covariance_m))
print("Param Messing", params_m, errors_m)
plt.plot(x_messing_fit * 1000, D_x_messing, 'rx', label="Messwerte")
plt.plot(x_linspace_fit_messing * 1000, D(x_linspace_fit_messing, *params_m),
         'b-', label="Regressionsgrade")
plt.xlabel(
    r"$L\cdot x^2- \frac{x^3}{3}$/$10^{-3}\,\si{\cubic\meter}$")
plt.ylabel(
    r"$D(x)$/$10^{-3}\,\si{\cubic\meter}$")
#plt.xlim(0, 53)
#plt.ylim(0, 5.6)
plt.axis('tight')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b.pdf')
ascii.write([x_messing * 100, d_messing - 0.5, dx_messing - 0.5, D_x_messing], 'Messdaten/messing_einseitig.tex', format="latex",
            names=["messpunkt x ", "D_0", "D_belastet", "D diff"])
a_messing = ufloat(params_m[0], errors_m[0])
E_messing = F_messing / (2 * a_messing * I_messing)
print("E Messing", E_messing)

#####################alu#################
x_linspace = np.linspace(0, 49) / 100

m_alu, d_alu, dx_alu = np.genfromtxt("Messdaten/a.txt", unpack=True)
D_x_alu = dx_alu - d_alu
x_alu = m_alu + 3
L_stab_alueingespannt = 0.54
mlast_alu = 767.4
L_stab_alu = 60
d_stab_alu = 0.01
v_stab_alu = d_stab_alu**2 * L_stab_alu
m_stab_alu = 167.1
pdichte_alu = m_stab_alu / v_stab_alu
print("Dichte Stab rechteckig", pdichte_alu)
dichte_lit_alu = 2.712  # in g/m^3
print("Dichte Alu Literatur", dichte_lit_alu)


x_alu = x_alu / 100
x_alu_fit = (L_stab_alueingespannt * x_alu**2 - ((x_alu**3) / 3))
x_linspace_fit_alu = (L_stab_alueingespannt *
                      x_linspace**2 - ((x_linspace**3) / 3))

ascii.write([x_alu * 100, d_alu - 1, dx_alu - 1, D_x_alu], 'Messdaten/alu_einseitig.tex', format="latex",
            names=["messpunkt x ", "D_0", "D_belastet", "D diff"])


plt.clf()
params, covariance = curve_fit(D, x_alu_fit, D_x_alu)
errors = np.sqrt(np.diag(covariance))
print("Param Alu", params, errors)
plt.plot(x_alu_fit * 1000, D_x_alu, 'rx', label="Messwerte")
plt.plot(x_linspace_fit_alu * 1000, D(x_linspace_fit_alu, *params),
         'b-', label="Regressionsgrade")
plt.xlabel(
    r"$L\cdot x^2- \frac{x^3}{3}$/$10^{-3}\,\si{\cubic\meter}$")
plt.ylabel(
    r"$D(x)$/$10^{-3}\,\si{\cubic\meter}$")
#plt.xlim(0, 53)
#plt.ylim(0, 5.75)
plt.axis('tight')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/a.pdf')
a_alu = ufloat(params[0], errors[0])
F_alu = mlast_alu * g
I_alu = d_stab_alu**4 / 12
print("I alu=", I_alu)
E_alu = F_alu / (2 * a_alu * I_alu)
print("E alu=", E_alu)
