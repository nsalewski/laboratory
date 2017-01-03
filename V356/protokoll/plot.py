import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
n,w=np.genfromtxt("Messdaten/c_1.txt",unpack=True)
n=n+1
theta=(n*np.pi)/14
kreisfrequenz=w*2*np.pi
phase=kreisfrequenz/theta
L=1.217*1/10**4
C=20.13*1/10**9
def theorie(f):
    return ((2*np.pi*f)/(np.arccos(1-((1/2)*L*C*(2*np.pi*f)**2))))

ascii.write([n, theta,w,w*2*np.pi,phase], 'Messdaten/tab_c.tex', format="latex",
            names=['n','theta','frequenz','kreis','Phase'])
plt.plot(w, phase, 'rx', label="Messwerte")
plt.plot(w, theorie(w), 'b-', label="Theoriekurve")

plt.ylabel(r'Phasengeschwindigkeit')
plt.xlabel(r'Frequenz')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/Phasengeschwindigkeit.pdf')
