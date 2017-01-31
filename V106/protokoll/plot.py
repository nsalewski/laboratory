import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

ascii.write([diff, diff * 2], 'Messdaten/entartung_wir.tex', format='latex')

plt.plot(v, g_noms, 'ro', label="Messwerte")
plt.plot(v, ausgleichsgrade(
    v, *paramsdeltav), 'b-', label="Regressionsgrade")
plt.xlabel(
    r"Geschwindigkeit des Reflektors $v_{\mathrm{R}}$/$\si{\meter\per\second}$")
plt.ylabel(
    r"Frequenzänderung $\Delta \nu$/$\si{Hz}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/schwebung.pdf')
