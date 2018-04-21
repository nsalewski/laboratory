#!usr/bin/env python
#coding:utf8
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import textable
import scipy.constants as const
import math as math
def kernspin(g_j,g_f):
    return(((g_j)/(2*g_f))-0.5)
def theorie(x,m,b):
    return m*x+b
def helmholtz(n,r,i):
    return((const.mu_0)*(8*i*n)/(r*np.sqrt(125)))
#Daten
rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)
vertical=0.231

#SI
rf_theo=1000*np.linspace(-10,1150)
rf=1000*rf
horizontal_1/=1000
horizontal_2/=1000
gj=2.0023
#B-Felder
B_vert=helmholtz(20,0.11735, vertical)
B_sweep_1=helmholtz(11,0.1639,peak_1)
B_sweep_2=helmholtz(11,0.1639,peak_2)
B_Horizontal_1=helmholtz(154, 0.1579,horizontal_1)
B_Horizontal_2=helmholtz(154,0.1579,horizontal_2)
B1=B_sweep_1+B_Horizontal_1
B2=B_sweep_2+B_Horizontal_2

#Tabellen
textable.latex_tab(data=[rf/1000,(horizontal_1*10**3).astype(int),(horizontal_2*10**3).astype(int),(peak_1*10**3).astype(int),(peak_2*10**3).astype(int)],names=[r"RF-Wechselfeld/$\si{\kilo\hertz}$",r"$I_{\mathrm{Horizontal,1}}/\si{\milli \ampere}$",r"$I_{\mathrm{Horizontal,2}}/\si{\milli \ampere}$",r"$I_{\mathrm{Sweep,1}}/\si{\milli \ampere}$",r"$I_{\mathrm{Sweep,2}}/\si{\milli \ampere}$"], filename=r"latex_tables/current.tex",caption=r"Im Experiment gemessene Ströme der Sweep-Spule und der Horizontalfeldspule für die Transparenzminima beider Isotope sowie die Frequenz des angelegten RF-Wechselfelds ",label=r"tab:current")

textable.latex_tab(data=[rf/1000,np.round(B_sweep_1*10**6,2),np.round(B_sweep_2*10**6,2),np.round(B_Horizontal_1*10**6,2),np.round(B_Horizontal_2*10**6,2)], names=[r"RF-Wechselfeld/$\si{\kilo\hertz}$",r"$B_{\mathrm{Horizontal,1}}/10^{-6}\si{\tesla}$",r"$B_{\mathrm{Horizontal,2}}/10^{-6}\si{\tesla}$",r"$B_{\mathrm{Sweep,1}}/10^{-6}\si{\tesla}$",r"$B_{\mathrm{Sweep,2}}/10^{-6}\si{\tesla}$"], filename=r"latex_tables/fields.tex",caption=r"Aus den gemessenen Strömen berechnete B-Felder für die Horizontalfeldspule und die Sweep-Spule in den Transparenzminima beider Isotope",label=r"tab:fields")


#Ausgleichsrechnung
params1, covariance1 = curve_fit(theorie,rf,B1)
errors1 = np.sqrt(np.diag(covariance1))

params2, covariance2 = curve_fit(theorie,rf,B2)
errors2 = np.sqrt(np.diag(covariance2))

a1=ufloat(params1[0],errors1[0])
a2=ufloat(params2[0],errors2[0])
erd1=ufloat(params1[1],errors1[1])
erd2=ufloat(params2[1],errors2[1])


erdges=np.mean([erd1,erd2])
g1=(4*np.pi*const.m_e)/(a1*const.e)
g2=(4*np.pi*const.m_e)/(a2*const.e)
print("I1: ",'{:.3f}'.format(kernspin(gj,g1)))
print("I2: ",'{:.3f}'.format(kernspin(gj,g2)))
#Ausgaben
print("B-Feld Vertikal=",np.round(B_vert*10**6,2))
print("Steigungsparam;Horizontales lokales Erdmagnetfeld1;g1",'{:.2f}'.format(a1*10**12),'{:.2f}'.format(erd1*10**6),'{:.3f}'.format(g1))
print("a2;Horizontales lokales Erdmagnetfeld2; g2",'{:.2f}'.format(a2*10**12),'{:.2f}'.format(erd2*10**6),'{:.3f}'.format(g2))
print("Horizontales lokales Erdmagnetfeld",'{:.1f}'.format(erdges*10**6))

#Amplitudenverhältnis
a1=ufloat(287,0.02*287)
a2=ufloat(496,0.02*496)
print(a1,a2,a1/a2)

#quadratischer Zeeman
mu_b=9.274*10**(-24)
def zeeman(gf,mf,ehy,b):
    return (((gf**2)*(mu_b**2)*(b**2)*((1-2*mf)/(ehy))))
print("Zeeman1",(g1*mu_b*250*10**(-6)),zeeman(g1,0,2.01*10**(-24),250*10**(-6)))
print("Zeeman2",(g2*mu_b*250*10**(-6)),zeeman(g2,0,4.53*10**(-24),250*10**(-6)))

#Plot
plt.plot(rf/1000,B1*10**6, 'ro', label="Messwerte 1. Minimum")
plt.plot(rf/1000,B2*10**6, 'rx', label="Messwerte 2. Minimum")
plt.plot(rf_theo/1000, theorie(rf_theo, *params1)*10**6, 'b-', label="Regressionsgrade 1. Minimum")
plt.plot(rf_theo/1000, theorie(rf_theo, *params2)*10**6, 'k-', label="Regressionsgrade 2.Minimum")
plt.xlim(0,1100)
plt.ylabel(r"Gesamtes Horizontales Magnetfeld $B_{\mathrm{ges}}/10^{-6}\si{\tesla}$")
plt.xlabel(r"Frequenz des RF-Felds $\nu/\si{\kilo\hertz}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('pictures/lin_regress.pdf')
