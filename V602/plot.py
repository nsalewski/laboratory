import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
Z=[29,29,29,30,32,35,37,38,40,41,79,79]
Z=np.asarray(Z)
E=[8.048,8.028,8.905,9.65,11.1031,13.4737,15.1997,16.1046,17.9976,18.9856,13.739,11.925]
E=np.asarray(E)
E=E*1000
print(E)
enull=1.6021766208*10**(-19)
R=13.6
r=10973731.568508
c=299792458
d=201.4*10**(-12)
h=ufloat(6.626070040*10**(-34),0.000000081*10**(-34))
Rzwo=r*h.nominal_value*c
R=Rzwo
E=E*enull
print(Z[0]-np.sqrt(E[0]/R))
theta=(360/(2*np.pi))*(np.arcsin(((c*h.nominal_value)/(2*d*E))))
sigma=Z-np.sqrt(E/R)
print(sigma)
print(theta)

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
