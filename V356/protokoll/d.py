import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
n1, u1 = np.genfromtxt("Messdaten/c_1.txt", unpack=True)
n2, u2 = np.genfromtxt("Messdaten/c_2.txt", unpack=True)
n3, u3 = np.genfromtxt("Messdaten/d.txt", unpack=True)
ascii.write([n1, u1, u2], 'Messdaten/tab_d.tex', format="latex",
            names=['Nummer des Kettenglieds', r'$U$/$\\si{\\volt}$ bei 1.Eigenschwingung', "$U$/$\\si{\\volt}$ bei 2.Eigenschwingung"])
plt.plot(n1, u1, 'ro--',
         label="Messdaten 1. Eigenschwingung (7.5 \,\si{\kilo\Hz})")
plt.ylabel(r'Spannung $U$/$\si{\volt}$')
plt.xlabel(r'Nr. des $LC$-Kettenglieds')
plt.xlim(-1, 14)
plt.ylim(0, 1.4)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/d.pdf')
plt.clf()

plt.plot(n2, u2, 'ro--',
         label='Messdaten 2. Eigenschwingung (14.7\,\si{\kilo\Hz})')
plt.ylabel(r'Spannung $U$/$\si{\volt}$')
plt.xlabel(r'Nr. des $LC$-Kettenglieds')
plt.legend(loc='best')
plt.xlim(-1, 14)
plt.tight_layout()
plt.savefig('Bilder/d2.pdf')
plt.clf()

ascii.write([n3, u3], 'Messdaten/tab_dlame.tex', format="latex",
            names=['Nummer des Kettenglieds', r'$U$/$\\si{\\volt}$ bei 1.Eigenschwingung und Abschlusswiderstand $R=Z$'])
plt.plot(n3, u3, 'ro--',
         label='Messdaten zur 1. Eigenschwingung (7.5\,\si{\kilo\Hz})')
plt.ylabel(r'Spannung $U$/$\si{\volt}$')
plt.xlabel(r'Nr. des $LC$-Kettenglieds')
plt.legend(loc='best')
plt.xlim(-1, 14)
plt.tight_layout()
plt.savefig('Bilder/d3.pdf')
