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

print('Mittelwert:', np.mean(Te), ' und Fehler des Mittelwerts: ',
      np.std(Te, ddof=1) / np.sqrt(len(Te)))


T_mittel = ufloat(np.mean(Te), np.std(Te, ddof=1) / np.sqrt(len(Te)))
Theta_ges = ufloat(1.3422, 0.0006) / 10000
magn_moment = ufloat(0.0339, 0.0004)
G = ufloat(90, 5) / 10**(9)
R = ufloat(86, 1.1) / 10**(6)
L = 58.5 / 100
B = ((((np.pi * 2) / T_mittel)**2 * Theta_ges) /
     (magn_moment)) - ((np.pi * G * R**4) / (2 * L * magn_moment))
print("Erdmagnetfeld= ", B)
