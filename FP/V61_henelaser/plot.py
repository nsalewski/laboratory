#!usr/bin/env python
#coding:utf8
#from __future__ import division
#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')
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
import matplotlib.ticker as tck
#Daten
#rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)
#Fehlmessungen:
#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])

# dec_points sets precision, i.e. dec_points[0]=2 will display 2 decimal places for all values in column 1
#Ausgleichsrechnung
#params1, covariance1 = curve_fit(theorie,rf,B1)
#errors1 = np.sqrt(np.diag(covariance1))

#vorbereitung
L=np.linspace(0,3)
r1=1.4
def f(r,L):
    return 1-L/r
plt.plot(L,f(r1,L)**2,label="beide Konkav")
plt.plot(L, f(r1,L), label="einer konkav, einer plan")
plt.plot([0,3],[1,1], label="obere Grenze")
plt.plot([0,3],[0,0], label="untere Grenze")
plt.legend(loc='best')
plt.ylabel(r"$g_1\cdot g_2$")
plt.xlabel(r"$L$/$\si{\meter}$")
plt.tight_layout()
plt.savefig("stabilitaet.pdf")
plt.clf()
#stabilität
def intens (phi,I_null,phi_null):
    return I_null*(np.sin(phi-phi_null))**2





phi,i=np.genfromtxt("daten/polarisation.txt",unpack=True)
phi=phi*2*np.pi/360


#Ausgleichsrechnung
params, covariance = curve_fit(intens,phi,i)
errors = np.sqrt(np.diag(covariance))
phi_theo=np.linspace(-20,160)
phi_theo=np.pi*2*phi_theo/360
f,ax=plt.subplots()
ax.plot(phi/np.pi,i,'bx',label="Messdaten")
ax.plot(phi_theo/np.pi, intens(phi_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I(\phi)$/$\si{\micro\ampere}$", xlabel=r"$\phi$/$\si{\radian}$")
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
plt.tight_layout()
plt.savefig("daten/polarisation.pdf")

def intens00(L, d_0, w,I_0):
    return I_0*np.exp(-2*((L-d_0)/w)**2)

vers, i00=np.genfromtxt("daten/tem_moden.txt",unpack=True)
params,covariance=curve_fit(intens00, vers, i00)
errors = np.sqrt(np.diag(covariance))
vers00_theo=np.linspace(0.9*np.min(vers), 1.1*np.max(vers))
f,ax=plt.subplots()
ax.plot(vers,i00,'bx',label="Messdaten")
ax.plot(vers00_theo, intens00(vers00_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{00}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/tem00.pdf")



def intens01(L, d_01, w1,I_01,d_02, w2,I_02):
    return I_01*np.exp(-2*((L-d_01)/w1)**2)+I_02*np.exp(-2*((L-d_02)/w2)**2)

vers01, i01=np.genfromtxt("daten/tem01.txt",unpack=True)
params,covariance=curve_fit(intens01, vers01, i01)
errors = np.sqrt(np.diag(covariance))
vers01_theo=np.linspace(0.9*np.min(vers01), 1.1*np.max(vers01))
f,ax=plt.subplots()
ax.plot(vers01,i01,'bx',label="Messdaten")
ax.plot(vers01_theo, intens01(vers01_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/tem01.pdf")


def intens01(L, d_01, w1,I_01,d_02, w2,I_02):
    return I_01*np.exp(-2*((L-d_01)/w1)**2)+I_02*np.exp(-2*((L-d_02)/w2)**2)

vers01, i01=np.genfromtxt("daten/tem01.txt",unpack=True)
params,covariance=curve_fit(intens01, vers01, i01)
errors = np.sqrt(np.diag(covariance))
vers01_theo=np.linspace(0.9*np.min(vers01), 1.1*np.max(vers01))
f,ax=plt.subplots()
ax.plot(vers01,i01,'bx',label="Messdaten")
ax.plot(vers01_theo, intens01(vers01_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/tem01.pdf")

def i_skk_fit(L,a,b,c):
    return a*L**2+b*L+c

def i_skp_fit(L,a,b,c):
    return a*L+b

vers_skk, i_skk=np.genfromtxt("daten/stabilitaet.txt",unpack=True)
params,covariance=curve_fit(i_skk_fit, vers_skk, i_skk)
errors = np.sqrt(np.diag(covariance))
vers_skk_theo=np.linspace(0.9*np.min(vers_skk), 1.1*np.max(vers_skk))
f,ax=plt.subplots()
ax.plot(vers_skk,i_skk,'bx',label="Messdaten")
ax.plot(vers_skk_theo, i_skk_fit(vers_skk_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/stabilitaetkk.pdf")

vers_skp, i_skp=np.genfromtxt("daten/stabilitaetpk.txt",unpack=True)
params,covariance=curve_fit(i_skp_fit, vers_skp, i_skp)
errors = np.sqrt(np.diag(covariance))
vers_skp_theo=np.linspace(0.9*np.min(vers_skp), 1.1*np.max(vers_skp))
f,ax=plt.subplots()
ax.plot(vers_skp,i_skp,'bx',label="Messdaten")
ax.plot(vers_skp_theo, i_skp_fit(vers_skp_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/stabilitaetkp.pdf")

def wl(a,n,d,L):
    return a/n*np.sin(np.arctan(d/L))
links,rechts,n=np.genfromtxt("daten/wellenlaenge.txt",unpack=True)
a=1/100*10**(-3)
L=28.8
for i, elem in enumerate(n):
    rechts[i]=wl(a,n[i],rechts[i],L)
    links[i]=wl(a,n[i],links[i],L)
print(rechts, links)
