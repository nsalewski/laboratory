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







#Daten
#rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)



#how to use textable

#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])

# dec_points sets precision, i.e. dec_points[0]=2 will display 2 decimal places for all values in column 1



#Ausgleichsrechnung
#params1, covariance1 = curve_fit(theorie,rf,B1)
#errors1 = np.sqrt(np.diag(covariance1))

#Plot
#plt.plot(rf/1000,B1*10**6, 'ro', label="Messwerte 1. Minimum")
#plt.plot(rf_theo/1000, theorie(rf_theo, *params1)*10**6, 'b-', label="Regressionsgrade 1. Minimum")
#plt.xlim(0,1100)
#plt.ylabel(r"$B_{\mathrm{ges}}/10^{-6}\si{\tesla}$")
#plt.xlabel(r"Frequenz des RF-Felds $\nu/\si{\kilo\hertz}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('pictures/lin_regress.pdf')
