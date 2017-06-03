import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
def form(t,c):
    return 0.5*c*t
def theorie(x,m,b):
    return m*x+b

s,t=np.genfromtxt("Messdaten/laufzeit.txt", unpack=True)
s=s*10**(-3)
t=t*10**(-6)
params, covariance = curve_fit(theorie,t,s)
errors = np.sqrt(np.diag(covariance))
c_halbe=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
c=2*c_halbe
print("c= ",c)
print("b= ",b)
plt.plot(t, s, 'rx', label="Messwerte")
t_l=np.linspace(20,90)
t_l=t_l*10**(-6)
plt.plot(t_l, theorie(t_l,*params), 'b-', label="Theoriekurve")

plt.ylabel(r"$s$/$\si{\milli\meter}")
plt.xlabel(r"$t$/$\si{\micro\second}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/ascan.pdf')
plt.clf()
s,t=np.genfromtxt("Messdaten/durchschall_messung_d.txt", unpack=True)
s=s*10**(-3)
t=t*10**(-6)
params, covariance = curve_fit(theorie,t,s)
errors = np.sqrt(np.diag(covariance))
c=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
print("c Durchschall= ",c)
print("b Durchschall= ",b)
plt.plot(t, s, 'rx', label="Messwerte")
t_ls=np.linspace(22,32)
t_ls=t*10**(-6)
plt.plot(t_ls, theorie(t_ls,*params), 'b-', label="Theoriekurve")
plt.ylabel(r"$s$/$\si{\milli\meter}")
plt.xlabel(r"$t$/$\si{\micro\second}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/d.pdf')
plt.clf()
s,a1,a2=np.genfromtxt("Messdaten/amplitude.txt", unpack=True)
s=s*10**(-3)

params, covariance = curve_fit(theorie,s,np.log(a1/a2))
errors = np.sqrt(np.diag(covariance))
a=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
c=2*c_halbe
print("a= ",-a)
print("b= ",b)
s_ls=np.linspace(25,88)
s_ls=s_ls*10**(-3)
plt.plot(s, np.log(a1/a2), 'rx', label="Messwerte")
plt.plot(s_ls, theorie(s_ls,*params), 'b-', label="Theoriekurve")

plt.ylabel(r"$\ln(\frac{I_0}{I(x)})$")
plt.xlabel(r"$s$/$\si{\milli\meter}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/amplitude.pdf')
plt.clf()

###############################################################################################
#Bestimmung der Abstände im Auge:
c_L = 2500
c_GK = 1410
n,t=np.genfromtxt('Messdaten/auge.txt', unpack=True)
t=t*10**(-6)

Strecke_Hornhaut_AnfangLinse = 1/2 * c_GK * (t[1]-t[0])
Strecke_Hornhaut_EndeLinse = 1/2 * c_L * (t[2]-t[1]) + Strecke_Hornhaut_AnfangLinse
Strecke_Hornhaut_Retina = 1/2 * c_GK * (t[3]-t[2]) + Strecke_Hornhaut_EndeLinse
print('Abstand Hornhaut Anfang Linse= ', Strecke_Hornhaut_AnfangLinse)
print('Avstand Hornhaut Ende Linse= ', Strecke_Hornhaut_EndeLinse)
print('Abstand Hornhaut Retina= ', Strecke_Hornhaut_Retina)
print(t)

###############################################################################################

#n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
#f=f*1000
#theta=(n*np.pi)/14
#w=f*2*np.pi
#L=1.217*1/10**3
#C=20.13*1/10**9
#thetaplot = np.linspace(0, 3)
#def theorie(theta):
#    return np.sqrt(2/(L*C)*(1-np.cos(theta)))
#ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
#            names=['n','frequenz','kreis','theta'])
#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#curve_fitting:
#def theorie(x,m,b):
#    return m*x+b
#ascii.write([np.sqrt(Ek),Z],'tab_007.tex',format='latex',names=["wurzel e","Z"])
#plt.plot(Z,np.sqrt(Ek), 'rx', label="Messwerte")
#params, covariance = curve_fit(theorie,Z,np.sqrt(Ek))
#errors = np.sqrt(np.diag(covariance))
#ryd=ufloat(params[0],errors[0])
