#!usr/bin/env python3
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
#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])
def manipulate(arr):
    for elem in range(len(arr)):
        if arr[elem-1]<180:
            arr[elem-1]=arr[elem-1]+180
        else: arr[elem-1]=arr[elem-1]-180
    return arr
def theorie(x,a,mu,b):
    return ((a*np.exp(-((x-mu)**2)/(b))))
def winkel(grad,sec):
    sec=sec*1/60
    grad=grad+sec
    return grad
def lin(x,a):
    return a*x
def eff_mass(a,B,N):
    return unp.sqrt(((e0)**3*N*B)/(8*np.pi**2*eps*c**3*n*a))
#daten importieren
b,z=np.genfromtxt("data/b_feld.txt",unpack=True)

f1,d1_hin,d1_hins,d1_rueck,d1_ruecks=np.genfromtxt("data/1_probe.txt",unpack=True)
f2,d2_hin,d2_hins,d2_rueck,d2_ruecks=np.genfromtxt("data/2_probe.txt",unpack=True)
f3,d3_hin,d3_hins,d3_rueck,d3_ruecks=np.genfromtxt("data/3_probe.txt",unpack=True)
f1=f1*10**(-6)
f2=f2*10**(-6)
f3=f3*10**(-6)
l1=1.296*10**(-3)
l2=1.36*10**(-3)
l3=5.11*10**(-3)

#bogensekunden addieren
grad1_hin=winkel(d1_hin,d1_hins)
grad1_rueck=winkel(d1_rueck,d1_ruecks)
grad2_hin=winkel(d2_hin,d2_hins)
grad2_rueck=winkel(d2_rueck,d2_ruecks)
grad3_hin=winkel(d3_hin,d3_hins)
grad3_rueck=winkel(d3_rueck,d3_ruecks)

#umrechnen auf gleichen Bezugspunkt
grad1_hin=manipulate(grad1_hin)
grad1_rueck=manipulate(grad1_rueck)
grad2_hin=manipulate(grad2_hin)
grad2_rueck=manipulate(grad2_rueck)
grad3_hin=manipulate(grad3_hin)
grad3_rueck=manipulate(grad3_rueck)

grad1=(1/(2*l1)*(grad1_rueck-grad1_hin)*2*np.pi/360)
grad2=(1/(2*l2)*(grad2_rueck-grad2_hin)*2*np.pi/360)
grad3=(1/(2*l3)*(grad3_rueck-grad3_hin)*2*np.pi/360)
#Berechnung delta theta
delta1=grad1-grad3
delta2=grad2-grad3
textable.latex_tab(data=[f1*10**6,grad3,grad1,grad2,delta1,delta2],names=[r"$\lambda$/$\si{\micro\meter}$",r"$\theta_{\mathrm{und}}$/$\si{\radian\per\meter}$",r"$\theta_{\mathrm{d1}}$/$\si{\radian\per\meter}$",r"$\theta_{\mathrm{d2}}$/$\si{\radian\per\meter}$",r"$\Delta \theta_{\mathrm{d1}}$/$\si{\radian\per\meter}$",r"$\Delta \theta_{\mathrm{d2}}$/$\si{\radian\per\meter}$"], filename=r"tables/eff_mass.tex",caption=r"Werte der $\Delta \theta$ zwischen undotiertem und dotiertem $\ce{GaAs}$ zur Bestimmung der effektiven Masse der Kristallelektronen",label=r"eff_mass",dec_points=[2,2,2,2,2,2],tableformat=4.2)



#Tabellen theta
textable.latex_tab(data=[f1*10**6,grad1_hin,grad1_rueck,grad1],names=[r"$\lambda$/$\si{\micro\meter}$",r"$\theta_1$/$\si{\degree}$",r"$\theta_2$/$\si{\degree}$",r"$\theta$/$\si{\radian\per\meter}$"], filename=r"tables/probe1.tex",caption=r"Messwerte der Faraday-Rotation für die dotierte Probe $\ce{GaAs}_{d1}$",label=r"probe1",dec_points=[2,2,2,2],tableformat=4.2)
textable.latex_tab(data=[f2*10**6,grad2_hin,grad2_rueck,grad2],names=[r"$\lambda$/$\si{\micro\meter}$",r"$\theta_1$/$\si{\degree}$",r"$\theta_2$/$\si{\degree}$",r"$\theta$/$\si{\radian\per\meter}$"], filename=r"tables/probe2.tex",caption=r"Messwerte der Faraday-Rotation für die dotierte Probe $\ce{GaAs}_{d2}$",label=r"probe2",dec_points=[2,2,2,2],tableformat=4.2)
textable.latex_tab(data=[f3*10**6,grad3_hin,grad3_rueck,grad3],names=[r"$\lambda$/$\si{\micro\meter}$",r"$\theta_1$/$\si{\degree}$",r"$\theta_2$/$\si{\degree}$",r"$\theta$/$\si{\radian\per\meter}$"], filename=r"tables/probe3.tex",caption=r"Messwerte der Faraday-Rotation für die undotierte Probe $\ce{GaAs}_{und}$",label=r"probe3",dec_points=[2,2,2,2],tableformat=4.2)
#Tabelle Magnetfeld
textable.latex_tab(data=[z-3.1,b],names=[r"$z$/$\si{\centi\meter}$",r"$B$/$\si{\milli\tesla}$"], filename=r"tables/magnetfeld.tex",caption=r"Messung des Magnetfelds in Abhängigkeit zum Ort $z$ (Probe ist etwa bei $\SI{3.1}{\centi\meter}$ platziert)",label=r"magnetfeld",dec_points=[2,0],tableformat=3.2)


z_theo=np.linspace(0,6,50)
#Ausgleichsrechnung Magnetfeld
params, covariance = curve_fit(theorie,z-3.1,b)
errors = np.sqrt(np.diag(covariance))
print(params,errors)
print("Erwartungswert",params[1],errors[1])
delta1_calc=np.delete(delta1,[0,3,7])
f1_calc1=np.delete(f1,[0,3,7])
delta2_calc=np.delete(delta2,[6,7])
f1_calc2=np.delete(f1,[6,7])


#lin regress delta

paramsd1, covarianced1 = curve_fit(lin,(f1_calc1**2),delta1_calc*10**(-6))
errorsd1 = np.sqrt(np.diag(covarianced1))
paramsd2, covarianced2 = curve_fit(lin,(f1_calc2)**2,delta2_calc*10**(-6))
errorsd2 = np.sqrt(np.diag(covarianced2))

a1=ufloat(paramsd1[0],errorsd1[0])*10**(6)
a2=ufloat(paramsd2[0],errorsd2[0])*10**(6)

n=3.3
e0=const.e
eps=const.epsilon_0
c=const.c
B=377.5*10**(-3)
print("Delta_1 Steigung", a1)
print("Delta_2 Steigung", a2)
print("Effektive Masse 1",eff_mass(a1,B,2.8*10**18*10**6),eff_mass(a1,B,2.8*10**18*10**6)/const.m_e)
print("Effektive Masse 2",eff_mass(a2,B,1.2*10**18*10**6),eff_mass(a2,B,1.2*10**18*10**6)/const.m_e)

#Plot Magnetfeld
plt.plot((params[1],params[1]),(-20,400), 'r--', label="Erwartungswert \n der Normalverteilung")
plt.plot(z-3.1,b, 'rx', label="Messwerte $B$")
plt.ylabel(r"$B/\si{\milli\tesla}$")
plt.xlabel(r"z/\si{\centi\meter}")
plt.legend(loc='best')
plt.ylim(-20,400)
axis.labels()
plt.tight_layout()
plt.savefig('pictures/B_feld.pdf')
plt.clf()

#Plot theta
plt.plot(f1*10**6,grad1, 'rx', label=r"Messwerte $\theta_{\mathrm{d1}}$")
plt.plot(f2*10**6,grad2, 'gx', label=r"Messwerte $\theta_{\mathrm{d2}}$")
plt.plot(f3*10**6,grad3, 'bx', label=r"Messwerte $\theta_{\mathrm{und}}$")
plt.ylabel(r"$\theta$/$\si{\radian\per\meter}")
plt.xlabel(r"$\lambda$/$\si{\micro\meter}$")
plt.legend(loc='lower right')
plt.tight_layout()
axis.labels()
plt.xlim(1,3.5)
plt.savefig('pictures/winkel_gg_wellenlaenge.pdf')
plt.clf()


f_theo=np.linspace(0,np.max(f1)+0.1*np.max(f1))
#plot delta
plt.plot((f1)**2*10**11,delta1, 'rx', label=r"$\Delta \theta_{\mathrm{d1}}$")
plt.plot((f_theo)**2*10**11,lin((f_theo)**2,*paramsd1*10**6), 'b-', label="Ausgleichsgrade")
plt.ylabel(r"$\Delta \theta_{\mathrm{d1}}$/$\si{\radian\per\meter}$")
plt.xlabel(r"$\lambda^{2}$/$\si{\square\meter}\cdot \num{e-11}$")
plt.legend(loc='best')
axis.labels()
plt.xlim(0,1.1)
plt.tight_layout()
plt.savefig('pictures/delta1.pdf')
plt.clf()
plt.plot((f1)**2*10**11,delta2, 'rx', label=r"$\Delta \theta_{\mathrm{d2}}$")
plt.plot((f_theo)**2*10**11,lin(f_theo**2,*paramsd2*10**6), 'b-', label="Ausgleichsgrade")
plt.ylabel(r"$\Delta \theta_{\mathrm{d2}}$/$\si{\radian\per\meter}$")
plt.xlabel(r"$\lambda^{2}$/$\si{\square\meter}\cdot\num{e-11}$")
axis.labels()
plt.legend(loc='best')
plt.tight_layout()
plt.xlim(0,1.1)

plt.savefig('pictures/delta2.pdf')
plt.clf()
