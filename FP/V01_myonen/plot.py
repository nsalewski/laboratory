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
from modules.plot import axislabel as axis
#Daten
#rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)
T=[1, 2, 3, 4, 5, 6, 7, 8, 9]
Kanal=[35, 57.04, 79.83, 102, 124, 146, 168, 190, 212, ]
texable.latex_tab(data=[T,2*Kanal], names[r'Kanalnummer', r'Doppelimpulsabstand T / $\si{\second}$'], filename=r'eichi.tex', caption=r"Kanalnummer in Abhängigkeit des Doppelimpulsabstandes $T$ für die Zeiteichung der Apparatur.", label=r"tab:eichi", dec_points=[2,0])

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
x = np.linspace(-2, 2, 100)
plt.plot(x, x**2, 'b-', label="test")
axis.labels()
plt.savefig('axislabel1.pdf')

x = np.linspace(-2, 2, 100)
plt.plot(x, x**2, 'b-', label="test")
plt.savefig('axislabel2.pdf')
x = np.linspace(-2, 2, 100)
plt.plot(x, x**2, 'b-', label="test")
plt.savefig('axislabel3.pdf')
