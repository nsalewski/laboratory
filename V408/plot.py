import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
f=f*1000
theta=(n*np.pi)/14
w=f*2*np.pi
L=1.217*1/10**3
C=20.13*1/10**9
thetaplot = np.linspace(0, 3)

def theorie(theta):
    return np.sqrt(2/(L*C)*(1-np.cos(theta)))

ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
            names=['n','frequenz','kreis','theta'])


plt.plot(theta, w/1000, 'rx', label="Messwerte")
plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")

plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
plt.xlabel(r"$\theta/\si{\radian}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b1.pdf')
