import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants

n, Te = np.genfromtxt("Messdaten/c.txt", unpack=True)
ascii.write([n, Te], 'Messdaten/c.tex', format="latex",
            names=["Anzahl", "Periode"])

print('Mittelwert:', np.mean(Te), ' und Fehler des Mittelwerts: ', np.std(Te, ddof=1)/np.sqrt(len(Te)))
