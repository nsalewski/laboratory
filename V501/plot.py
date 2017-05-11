import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

##########################################################################################
# E-Feld
n_, v_säge = np.genfromtxt("Messdaten/frequenzsaege.txt",unpack=True)
ascii.write([n_, v_säge], 'Messdaten/tab_saegi.tex', format="latex", names=['Frequenzverhältnis','frequenz'])

D, Ud400, Ud300, Ud200 = np.genfromtxt("Messdaten/efeld.txt",unpack=True)
ascii.write([D*2.54, Ud400, Ud300, Ud200], 'Messdaten/tab_efeld.tex', format="latex")



#########################################################################################
# B-Feld

I250, D_, I450  = np.genfromtxt("Messdaten/messdaten502a.txt",unpack=True)
ascii.write([D_*2.54, I250, I450], 'Messdaten/tab_bfeld.tex', format="latex")

def y(x, m, b):
    return m * x + b

params, covariance = curve_fit(y, 4*np.pi*10**(-7)*8/np.sqrt(125)*20*I250/0.282, D_/(D_**2+0.143**2))
errors = np.sqrt(np.diag(covariance))

print('m = ', params[0], '+/-', errors[0])
print('b = ', params[1], '+/-', errors[1])

plt.plot(np.linspace(0,0.0002), params[0]*np.linspace(0,0.0002)+params[1], 'b-',label='Ausgleichsgerade')
plt.plot(4*np.pi*10**(-7)*8/np.sqrt(125)*20*I250/0.282,D_/(D_**2+0.143**2) , 'rx', label='Messwerte')
plt.ylabel(r"$\frac{D}{D^2+L^2}$/$\si{\per\meter}$")
plt.xlabel(r"$B$/$\si{\tesla}$")
plt.tight_layout()
plt.savefig('Messdaten/plotbfeld.pdf')
















#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#
