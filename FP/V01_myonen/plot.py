#!usr/bin/env python
#coding:utf8
from __future__ import division
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
from modules.table import textable
import scipy.constants as const
import math as math
from modules.plot import axislabel as axis
def lin(x,a,b):
    return a*x+b
#Daten
#rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)
#Fehlmessungen:
n_start=1998242
zeit_1=94990
zeit_2=94991
zeit_ges=ufloat(np.mean([zeit_1,zeit_2]), stats.sem([zeit_1,zeit_2]))
f=n_start/zeit_ges
Ts=20*10**(-6)
fehl=f*Ts*n_start*unp.exp(f*Ts)
print(fehl, zeit_ges)
kanalzahl=442
u=fehl/kanalzahl
print("Frequenz", f, "Untergrundrate",u)
#Plateau:
delay,counts=np.genfromtxt("Daten/plateau.txt", unpack=True)
textable.latex_tab(data=[delay,counts], names=[r'Verzögerung $t_2$/$\si{\nano\second}$', r'Counts $N(t)$/$\SI{60}{\second}$'], filename=r'plateau.tex', caption=r"Messung zur Bestimmung der Verzögerungszeit der Apparatur und der Bestimmung der Koinzidenzbreite der Koinzidenz.", label=r"plateau", dec_points=[1,0],tableformat=3.3)
#fit links

params_lin, covariance_lin = curve_fit(lin,delay[:7],counts[:7])
errors_lin = np.sqrt(np.diag(covariance_lin))
#fit rechts
params_r, covariance_r = curve_fit(lin,delay[15:],counts[15:])
errors_r = np.sqrt(np.diag(covariance_r))
print("params lin", params_lin, errors_lin)
print("params r", params_r, errors_r)
plateau=np.mean(counts[7:15])
print("Plateau",plateau)
x1_halb=(plateau/2-params_lin[1])/params_lin[0]
x2_halb=(plateau/2-params_r[1])/params_r[0]
print("koinzidenzbreite", x2_halb-x1_halb,x1_halb,x2_halb)

plt.plot(delay,counts, 'rx', label="Daten")
plt.errorbar(delay, counts, yerr=np.sqrt(counts),fmt='none')
plt.plot([delay[7], delay[14]], [plateau,plateau], label="Plateau")
plt.plot([x1_halb, x2_halb], [plateau/2,plateau/2], label="Halbwertsbreite")
plt.plot(delay[:7], lin(delay[:7],*params_lin), label="Regression links")
plt.plot(delay[15:], lin(delay[15:],*params_r), label="Regression rechts")

plt.ylabel(r" Zählrate $N(t)$ / $\SI{60}{\second}$")
plt.xlabel(r"Verzögerung $\Delta t$/\si{\nano\second}")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Daten/plateau.pdf')
plt.clf()



T=[1, 2, 3, 4, 5, 6, 7, 8, 9]
T=np.array(T)
T=T*10**(-6)
Kanal=[23,45*3417/(3417+143)+46*143/(3417+143),67*710/(3546+710)+68*3546/(3546+710), 90, 112, 134, 156, 178, 200]
anzahl=[4432,3417+143,3546+710,4183,3657,5148,3805,3259,26001]
#print(2*Kanal) ist liste. Unterschied zu np.array.
Kanal=np.array(Kanal)
textable.latex_tab(data=[Kanal,T*10**6,anzahl], names=[r'Kanalnummer', r'Doppelimpulsabstand T / $\si{\micro\second}$',r'Anzahl an gemessenen Impulsen'], filename=r'eichi.tex', caption=r"Kanalnummer in Abhängigkeit des Doppelimpulsabstandes $T$ für die Zeiteichung der Apparatur.", label=r"eichi", dec_points=[2,0,0],tableformat=3.3)
#textable war falsch geschrieben und hinter names fehlte ein "=". Das tableformat was ich noch eingefügt habe, ist optional, die erste zahl steuert, wieviel platz global vor dem komma in der tabelle reserviert wird und die zweite, wieviel danach. für lange zahlen manchmal sehr sinnvoll.
#how to use textable

#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])

# dec_points sets precision, i.e. dec_points[0]=2 will display 2 decimal places for all values in column 1
def f(x, m, b):
    return m*x+b
params, covariance = curve_fit(f,Kanal,T)
errors = np.sqrt(np.diag(covariance))
print('m = ', params[0], 'pm', errors[0])
print('b = ', params[1], 'pm', errors[1])
Kanal_regress=np.linspace(0,220)
plt.plot(Kanal,T*10**6, 'rx', label="Daten")
plt.plot(Kanal_regress, f(Kanal_regress, *params)*10**6, 'b-', label="Regressionsgrade")
plt.ylabel(r"Impulsabstand $t$ / $\si{\micro\second}$")
plt.xlabel(r"Kanalnummer")
plt.xlim(0,220)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Daten/kalibrate.pdf')
plt.clf()
nmbrs = np.genfromtxt("Daten/messung.txt",unpack=True)
print('*******************')
print(sum(nmbrs))
print('*******************')

kanal = np.linspace(3, 444, 442)
t = kanal * params[0]+params[1]
def N(t, NO, l,b):
    return NO * np.exp(-l*t) + b
nmbrs_fit=[]
nmbrs_fit=np.array(nmbrs_fit)
nmbrs_fit=np.append(nmbrs_fit,nmbrs[2])
nmbrs_fit=np.append(nmbrs_fit,nmbrs[4:6])
nmbrs_fit=np.append(nmbrs_fit,nmbrs[7])
nmbrs_fit=np.append(nmbrs_fit,nmbrs[9:11])
nmbrs_fit=np.append(nmbrs_fit,nmbrs[14:])
nmbrs_ignore=[]
nmbrs_ignore=np.array(nmbrs_ignore)
nmbrs_ignore=np.append(nmbrs_ignore,nmbrs[:2])
nmbrs_ignore=np.append(nmbrs_ignore,nmbrs[3])
nmbrs_ignore=np.append(nmbrs_ignore,nmbrs[6])
nmbrs_ignore=np.append(nmbrs_ignore,nmbrs[8])

nmbrs_ignore=np.append(nmbrs_ignore,nmbrs[11:14])
t_ignore=[]
t_ignore=np.array(t_ignore)
t_ignore=np.append(t_ignore,t[:2])
t_ignore=np.append(t_ignore,t[3])
t_ignore=np.append(t_ignore,t[6])
t_ignore=np.append(t_ignore,t[8])

t_ignore=np.append(t_ignore,t[11:14])

t_fit=[]
t_fit=np.array(t_fit)
t_fit=np.append(t_fit,t[2])
t_fit=np.append(t_fit,t[4:6])
t_fit=np.append(t_fit,t[7])
t_fit=np.append(t_fit,t[9:11])
t_fit=np.append(t_fit,t[14:])
print(nmbrs_fit)

params1, covariance1 = curve_fit(N,t_fit,nmbrs_fit)

errors1 = np.sqrt(np.diag(covariance1))
print('NO = ', params1[0], 'pm', errors1[0])
print('l = ', params1[1], 'pm', errors1[1])
print('b = ', params1[2], 'pm', errors1[2])

arr = ufloat(params1[1], errors1[1])
print(1/arr)

#Ausgleichsrechnung
#params1, covariance1 = curve_fit(theorie,rf,B1)
#errors1 = np.sqrt(np.diag(covariance1))

t_plot = np.linspace(0, 21)*10**(-6)
#Plot
plt.plot(t_fit*10**6,nmbrs_fit, 'r+', label="Daten")
plt.plot(t_ignore*10**6, nmbrs_ignore, 'kx',label="Nicht betrachtete Daten")
plt.plot(t_plot*10**6, N(t_plot, *params1), 'b--', label="Theoriekurve")
#plt.xlim(0,1100)
plt.ylabel(r"$N(t)$")
plt.xlabel(r"$t$ / $\si{\micro\second}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Daten/regress.pdf')
