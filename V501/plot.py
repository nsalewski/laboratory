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
















#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#
