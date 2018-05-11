#!usr/bin/env python
#coding:utf8
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
from modules.table import textable
import scipy.constants as const
import math as math

#Konstanten

h = 6.626070040 *10**(-34) #Plancksche Wirkungsquantum
c = 299792458#Lichtgeschw
mu_B = 927.4009994 * 10**(-26) #Bohrsches Magneton
rot_lambda_D = 48.913 * 10**(-12)
blau_lambda_D = 26.952 * 10**(-12)

##############################################################################################
# Eichen
I, B = np.genfromtxt('data/magnetfeld.txt', unpack=True)

def polynom(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

params, covariance = curve_fit(polynom, I, B)
errors = np.sqrt(np.diag(covariance))
interval_eich = np.linspace(0, 16)

plt.plot(I, B, 'rx', label='Messwerte')
plt.plot(interval_eich, polynom(interval_eich, *params), 'b-', label='Ausgleichspolynom')
plt.xlabel(r'Feldstrom $I / \si{\ampere}$')
plt.ylabel(r'Magnetfeld $B / \si{\milli\tesla}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('pictures/eichi.pdf')
plt.clf()

##############################################################################################

lambd = 643.8*10**(-9)
B_rot = 652*10**(-3)
Ds = [148, 136, 130, 132, 132, 126, 120, 122, 122, 120]
ds = [72, 72, 68, 66, 70, 70, 66, 62, 62, 66]

D_s_mittel = ufloat(np.mean(Ds),np.std(Ds,ddof=1)/np.sqrt(len(Ds)))
d_s_mittel = ufloat(np.mean(ds),np.std(ds,ddof=1)/np.sqrt(len(ds)))
print('D_s Mittelwert: ', D_s_mittel)
print('d_s Mittelwert: ', d_s_mittel)
delta_lambda = 0.5 * d_s_mittel/D_s_mittel * rot_lambda_D
print('delta_lambda: ', delta_lambda)
g = (h*c* delta_lambda) / (lambd**2*mu_B*B_rot)
print('g: ', g)

###########################################################################

lambd = 480*10**(-9)
B_blau_sigma = 333*10**(-3)
blau_lambda_D = 26.952 * 10**(-12)
Ds = [126, 123, 117, 117, 114, 117, 111, 111, 112, 114 ]
ds = [60, 57, 61, 57, 60, 57, 54, 51, 57, 55]

D_s_mittel = ufloat(np.mean(Ds),np.std(Ds,ddof=1)/np.sqrt(len(Ds)))
d_s_mittel = ufloat(np.mean(ds),np.std(ds,ddof=1)/np.sqrt(len(ds)))
print('BLAU SIGMA :::: D_s Mittelwert: ', D_s_mittel)
print('d_s Mittelwert: ', d_s_mittel)
delta_lambda = 0.5 * d_s_mittel/D_s_mittel * blau_lambda_D
print('delta_lambda: ', delta_lambda)
g = (h*c* delta_lambda) / (lambd**2*mu_B*B_blau_sigma)
print('g: ', g)

########################################################################

lambd = 480*10**(-9)
B_blau_pi = 1001*10**(-3)
blau_lambda_D = 26.952 * 10**(-12)
Ds = [126, 120, 111, 117, 111, 105, 108, 102, 108, 100 ]
ds = [48, 51, 45, 46, 45, 48, 36, 42, 36, 35]

D_s_mittel = ufloat(np.mean(Ds),np.std(Ds,ddof=1)/np.sqrt(len(Ds)))
d_s_mittel = ufloat(np.mean(ds),np.std(ds,ddof=1)/np.sqrt(len(ds)))
print('BLAU PI :::: D_s Mittelwert: ', D_s_mittel)
print('d_s Mittelwert: ', d_s_mittel)
delta_lambda = 0.5 * d_s_mittel/D_s_mittel * blau_lambda_D
print('delta_lambda: ', delta_lambda)
g = (h*c* delta_lambda) / (lambd**2*mu_B*B_blau_pi)
print('g: ', g)


