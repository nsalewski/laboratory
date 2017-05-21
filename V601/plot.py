import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
n=np.genfromtxt("Messdaten/datenraumtemperatur.txt",unpack=True)
print(n)
c=np.mean(n)
m=np.std(n) / np.sqrt(len(n))
v1=ufloat(c,m)
v1=1/v1
print('Skala Raumtemperatur',v1)
p=np.genfromtxt("Messdaten/skalafranckhertz.txt",unpack=True)
c=np.mean(p)
m_2=np.std(p) / np.sqrt(len(p))
vhertz=ufloat(c,m_2)
print('Skala Franckhertz',vhertz)

vhertz=5/vhertz
print('Skala Franckhertz',vhertz)
e_theo = unp.uarray(1.6021766208*10**(-19), 0.0000000098*10**(-19))
abstand=np.genfromtxt("Messdaten/maximafranckhertz.txt", unpack=True)
ascii.write([abstand,abstand*vhertz.nominal_value],'Messdaten/franck.tex',format='latex')
c=np.mean(abstand*vhertz.nominal_value)
m_3=np.std(abstand*vhertz.nominal_value)/np.sqrt(len(abstand*vhertz.nominal_value))
print(c,m_3)
delta=ufloat(c,m_3)
e_theo = unp.uarray(1.6021766208*10**(-19), 0.0000000098*10**(-19))
delta=delta*e_theo
h=ufloat(6.626070040*10**(-34),0.000000081*10**(-34))
licht=299792458
lambdi=licht*h/delta
print(lambdi)
#n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
#f=f*1000
#theta=(n*np.pi)/14
#w=f*2*np.pi
#L=1.217*1/10**3
#C=20.13*1/10**9
#thetaplot = np.linspace(0, 3)
#
#def theorie(theta):
#    return np.sqrt(2/(L*C)*(1-np.cos(theta)))
#
#ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
#            names=['n','frequenz','kreis','theta'])
#
#
#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#
