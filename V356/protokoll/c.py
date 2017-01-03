import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
f=f*1000
ft=np.linspace(0,64.4)
ft=ft*1000
theta=(n*np.pi)/14
w=f*2*np.pi
phase=w/theta
L=1.217*1/10**3
C=20.13*1/10**9
def theorie(f):
    return ((2*np.pi*f)/(np.arccos(1-((1/2)*L*C*(2*np.pi*f)**2))))

ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2),np.round(phase/1000,2)], 'Messdaten/tab_c.tex', format="latex",
            names=['n','frequenz','kreis','theta','Phase'])


#divided by 1000 for proper formatting
plt.plot(f/1000, phase/1000, 'rx', label="Messwerte")
plt.plot(ft/1000, theorie(ft)/1000, 'b-', label="Theoriekurve")

plt.ylabel(r"$v_{\mathrm{Ph}}$/$\si{\kilo\metre\per\second}$")
plt.xlabel(r"$\omega/\si{\kilo\Hz}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/c.pdf')
