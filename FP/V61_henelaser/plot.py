#!usr/bin/env python
#coding:utf8
#from __future__ import division
#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')
from __future__ import division
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

#polarisation
def intens (phi,I_null,phi_null):
    return I_null*(np.sin(phi-phi_null))**2
phi,i=np.genfromtxt("daten/polarisation.txt",unpack=True)
phi=phi+20
phi=unp.uarray(phi*2*np.pi/360, len(phi)*[1*2*np.pi/360])
i=unp.uarray(i, len(i)*[0.1])
params, covariance = curve_fit(intens,unp.nominal_values(phi),unp.nominal_values(i))
errors = np.sqrt(np.diag(covariance))
print('****************************')
print("Polarisation")
print('I_0=%f(%f)'%(params[0],errors[0]),'\phi_0=%f(%f)'%(params[1], errors[1]))
print('****************************')

phi_theo=np.linspace(0,180)
phi_theo=np.pi*2*phi_theo/360
f,ax=plt.subplots()
ax.errorbar(unp.nominal_values(phi)/np.pi, unp.nominal_values(i),xerr=unp.std_devs(phi)/np.pi,yerr=unp.std_devs(i),fmt='bx',label="Messdaten")
ax.plot(phi_theo/np.pi, intens(phi_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I(\phi)$/$\si{\micro\ampere}$", xlabel=r"$\phi$/$\si{\radian}$")
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
plt.tight_layout()
plt.savefig("daten/polarisation.pdf")
textable.latex_tab(data=[phi,i],names=[r"Polarisationsmessung $\phi$/$\si{\degree}$",r"Intensität $I_{\mathrm{Pol}}$/$\si{\micro\ampere}$"], filename=r"daten/polarisation.tex",caption=r"Messdaten der Polarisationsmessung",label=r"polarisation",dec_points=[0,2])


#tem00
def intens00(L, d_0, w,I_0):
    return I_0*np.exp(-2*((L-d_0)/w)**2)

vers00, i00=np.genfromtxt("daten/tem_moden.txt",unpack=True)
vers00=unp.uarray(vers00-vers00[0],len(vers00)*[0.05])
i00=unp.uarray(i00, len(i00)*[0.2])
params,covariance=curve_fit(intens00, unp.nominal_values(vers00), unp.nominal_values(i00))
errors = np.sqrt(np.diag(covariance))
print('****************************')
print("TEM00")
print('d_0=%f(%f)'%(params[0],errors[0]),'w=%f(%f)'%(params[1], errors[1]),'I_0=%f(%f)'%(params[2], errors[2]))
print('****************************')
vers00_theo=np.linspace(0.9*np.min(unp.nominal_values(vers00)), 1.1*np.max(unp.nominal_values(vers00)))
f,ax=plt.subplots()
ax.errorbar(unp.nominal_values(vers00),unp.nominal_values(i00),xerr=unp.std_devs(vers00), yerr=unp.std_devs(i00),fmt='bx',label="Messdaten")
ax.plot(vers00_theo, intens00(vers00_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{00}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/tem00.pdf")
textable.latex_tab(data=[vers00,i00],names=[r"Verschiebung $\Delta L$/$\si{\centi\meter}$",r"Intensität $I_{\mathrm{00}}$/$\si{\micro\ampere}$"], filename=r"daten/i00.tex",caption=r"Messdaten der Messung der Grundmode $T_{00}$",label=r"t00",dec_points=[0,2],tableformat=3.3)


#tem01
def intens01(L, d_01, w1,I_01,d_02, w2,I_02):
    return I_01*np.exp(-2*((L-d_01)/w1)**2)+I_02*np.exp(-2*((L-d_02)/w2)**2)

vers01, i01=np.genfromtxt("daten/tem01.txt",unpack=True)
vers01=unp.uarray(vers01-vers01[0],len(vers01)*[0.05])
i01=unp.uarray(i01, len(i01)*[0.05])
params,covariance=curve_fit(intens01, unp.nominal_values(vers01), unp.nominal_values(i01),p0=[0.5,1.5,0.5,1.4,0.5,1.0])
errors = np.sqrt(np.diag(covariance))
print('****************************')
print("TEM01")
print('d_01=%f(%f)'%(params[0],errors[0]),'w1=%f(%f)'%(params[1], errors[1]),'I_01=%f(%f)'%(params[2], errors[2]),'d_02=%f(%f)'%(params[3],errors[3]),'w2=%f(%f)'%(params[4], errors[4]),'I_02=%f(%f)'%(params[5], errors[5]))
print('****************************')
vers01_theo=np.linspace(0.9*np.min(unp.nominal_values(vers01)), 1.1*np.max(unp.nominal_values(vers01)))
f,ax=plt.subplots()
ax.errorbar(unp.nominal_values(vers01),unp.nominal_values(i01),xerr=unp.std_devs(vers01), yerr=unp.std_devs(i01),fmt='bx',label="Messdaten")
ax.plot(vers01_theo, intens01(vers01_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\milli\meter}$")
plt.tight_layout()
plt.savefig("daten/tem01.pdf")
textable.latex_tab(data=[vers01,i01],names=[r"Verschiebung $\Delta L$/$\si{\centi\meter}$",r"Intensität $I_{\mathrm{01}}$/$\si{\micro\ampere}$"], filename=r"daten/i01.tex",caption=r"Messdaten der Messung der $T_{01}$-Mode",label=r"t01",dec_points=[0,2],tableformat=3.3)

#stabilität
def i_skk_fit(L,a,b,c):
    return a*L**2+b*L+c

def i_skp_fit(L,a,b):
    return a*L+b

vers_skk, i_skk=np.genfromtxt("daten/stabilitaet.txt",unpack=True)
vers_skk=unp.uarray(vers_skk+1,len(vers_skk)*[0.3])
i_skk=unp.uarray(i_skk, len(i_skk)*[0.2])
params,covariance=curve_fit(i_skk_fit, unp.nominal_values(vers_skk), unp.nominal_values(i_skk))
errors = np.sqrt(np.diag(covariance))
print('****************************')
print("Stabilitaet konkav-konkav: ")
print(params,errors)
print('a=%f(%f)'%(params[0],errors[0]),'b=%f(%f)'%(params[1], errors[1]),'c=%f(%f)'%(params[2], errors[2]))
print('****************************')
vers_skk_theo=np.linspace(0.9*np.min(unp.nominal_values(vers_skk)), 1.1*np.max(unp.nominal_values(vers_skk)))
f,ax=plt.subplots()
ax.errorbar(unp.nominal_values(vers_skk),unp.nominal_values(i_skk),xerr=unp.std_devs(vers_skk), yerr=unp.std_devs(i_skk),fmt='bx',label="Messdaten")
ax.plot(vers_skk_theo, i_skk_fit(vers_skk_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\centi\meter}$")
plt.tight_layout()
plt.savefig("daten/stabilitaetkk.pdf")
textable.latex_tab(data=[vers_skk,i_skk],names=[r"Resonatorlänge $ L$/$\si{\centi\meter}$",r"Intensität $I_{\mathrm{konkav-konkav}}$/$\si{\micro\ampere}$"], filename=r"daten/stabilitaetkk.tex",caption=r"Messdaten der Stabilitätsmessung für zwei konkave Resonatorspiegel",label=r"stabilitaetkk",dec_points=[0,2],tableformat=3.3)
vers_skp, i_skp=np.genfromtxt("daten/stabilitaetpk.txt",unpack=True)
vers_skp=unp.uarray(vers_skp+1,len(vers_skp)*[0.3])
i_skp=unp.uarray(i_skp, len(i_skp)*[0.025])
params,covariance=curve_fit(i_skp_fit, unp.nominal_values(vers_skp), unp.nominal_values(i_skp))
errors = np.sqrt(np.diag(covariance))
print('****************************')
print("Stabilitaet planar-konkav: ")
print('a=%f(%f)'%(params[0],errors[0]),'b=%f(%f)'%(params[1], errors[1]))
print('****************************')
vers_skp_theo=np.linspace(0.9*np.min(unp.nominal_values(vers_skp)), 1.1*np.max(unp.nominal_values(vers_skp)))
f,ax=plt.subplots()
ax.errorbar(unp.nominal_values(vers_skp),unp.nominal_values(i_skp),xerr=unp.std_devs(vers_skp),yerr=unp.std_devs(i_skp),fmt='bx',label="Messdaten")
ax.plot(vers_skp_theo, i_skp_fit(vers_skp_theo,*params),'r-', label="Fit")
ax.legend(loc='best')
ax.set(ylabel=r"$I_{01}(L)$/$\si{\micro\ampere}$", xlabel=r"$L$/$\si{\centi\meter}$")
plt.tight_layout()
plt.savefig("daten/stabilitaetkp.pdf")
textable.latex_tab(data=[vers_skp,i_skp],names=[r"Resonatorlänge $\Delta L$/$\si{\centi\meter}$",r"Intensität $I_{\mathrm{planar-konkav}}$/$\si{\micro\ampere}$"], filename=r"daten/stabilitaetpk.tex",caption=r"Messdaten der Stabilitätsmessung für einen konkaven und einen planaren Resonatorspiegel",label=r"stabilitaetpk",dec_points=[0,2],tableformat=3.3)
#wellenlaenge
def wl(a,n,d,L):
    return a/n*unp.sin(unp.arctan(d/L))*10**9
links,rechts,n=np.genfromtxt("daten/wellenlaenge.txt",unpack=True)
a=1/100*10**(-3)
L=ufloat(30.8,0.1)
err=[0.05,0.05,0.05]
links=unp.uarray(links, err)
rechts=unp.uarray(rechts, err)
rechtswl=unp.uarray(unp.nominal_values(wl(a,n,rechts,L)), unp.std_devs(wl(a,n,rechts,L)))
linkswl=unp.uarray(unp.nominal_values(wl(a,n,links,L)), unp.std_devs(wl(a,n,links,L)))
gesamt=unp.uarray([unp.nominal_values(rechtswl),unp.nominal_values(linkswl)],[unp.std_devs(rechtswl),unp.std_devs(linkswl)])
print('****************************')
print("Mittlere Wellenlänge: ",gesamt.mean())
print('****************************')
textable.latex_tab(data=[n,rechts,rechtswl,links,linkswl],names=[r"Nr. des Nebenmaxima",r"Abstand $r_{\mathrm{rechts}}$/$\si{\centi\meter}$",r"$\lambda_{\mathrm{rechts}}$/$\si{\nano\meter}$",r"Abstand $r_{\mathrm{links}}$/$\si{\centi\meter}$",r"$\lambda_{\mathrm{links}}$/$\si{\nano\meter}$"], filename=r"daten/wellenlaenge.tex",caption=r"Gemessene Abstände $r_i$ der Nebenmaxima zum Hauptmaximum und daraus berechnete Wellenlängen $\lambda_i$.",label=r"wellenlaenge",dec_points=[0,2,0,2,0], tableformat=3.2)
