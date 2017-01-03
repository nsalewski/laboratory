import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

n,f=np.genfromtxt("Messdaten/b_1.txt",unpack=True)
f=f*1000
theta=(n*np.pi)/14
w=f*2*np.pi
L=1.217*1/10**3
C1=20.13*1/10**9
C2=9.41*1/10**9
thetaplot = np.linspace(0, 2)

def theoriep(theta):
    return np.sqrt(1/L*(1/C1+1/C2)+1/L*np.sqrt((1/C1+1/C2)**2-(2*np.sin(theta))**2/(C1*C2)))

def theoriem(theta):
    return np.sqrt(1/L*(1/C1+1/C2)-1/L*np.sqrt((1/C1+1/C2)**2-(2*np.sin(theta))**2/(C1*C2)))


ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
            names=['n','frequenz','kreis','theta'])


plt.plot(theta, w/1000, 'rx', label="Messwerte")
plt.plot(thetaplot, theoriep(thetaplot)/1000, 'g-', label="Theoriekurve $\omega_+$")
plt.plot(thetaplot, theoriem(thetaplot)/1000, 'b-', label="Theoriekurve $\omega_-$")
plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
plt.xlabel(r"$\theta/\si{\radian}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b2.pdf')
