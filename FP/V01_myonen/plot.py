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
Kanal=[35, 57.04, 79.83, 102, 124, 146, 168, 190, 212]
#print(2*Kanal) ist liste. Unterschied zu np.array.
Kanal=np.array(Kanal)
textable.latex_tab(data=[T,Kanal*2], names=[r'Kanalnummer', r'Doppelimpulsabstand T / $\si{\second}$'], filename=r'eichi.tex', caption=r"Kanalnummer in Abhängigkeit des Doppelimpulsabstandes $T$ für die Zeiteichung der Apparatur.", label=r"tab:eichi", dec_points=[2,0],tableformat=3.3)
#textable war falsch geschrieben und hinter names fehlte ein "=". Das tableformat was ich noch eingefügt habe, ist optional, die erste zahl steuert, wieviel platz global vor dem komma in der tabelle reserviert wird und die zweite, wieviel danach. für lange zahlen manchmal sehr sinnvoll.
#how to use textable

#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])

# dec_points sets precision, i.e. dec_points[0]=2 will display 2 decimal places for all values in column 1
def f(x, m, b):
    return m*x+b
params, covariance = curve_fit(f,2*Kanal,T)
errors = np.sqrt(np.diag(covariance))
print('m = ', params[0], 'pm', errors[0])
print('b = ', params[1], 'pm', errors[1])

nmbrs = np.genfromtxt("Daten/messung.txt",unpack=True)
kanal = np.linspace(3, 444, 442)
t = kanal * params[0]
def N(t, NO, l):
    return NO * np.exp(-l*t) + 0.95


params1, covariance1 = curve_fit(N,t,nmbrs)
errors1 = np.sqrt(np.diag(covariance1))
print('NO = ', params1[0], 'pm', errors1[0])
print('l = ', params1[1], 'pm', errors1[1])
arr = ufloat(params1[1], errors1[1])
print(1/arr)

#Ausgleichsrechnung
#params1, covariance1 = curve_fit(theorie,rf,B1)
#errors1 = np.sqrt(np.diag(covariance1))

t_plot = np.linspace(0, 11)
#Plot
plt.plot(t,nmbrs, 'rx', label="Messwerte")
plt.plot(t_plot, N(t_plot, *params1), 'b-', label="Theoriekurve")
#plt.xlim(0,1100)
plt.ylabel(r"$N(t)$")
plt.xlabel(r"$t$ / $\si{\micro\second}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Daten/regress.pdf')

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
