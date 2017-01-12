import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
m_k = ufloat(512.2, 0.2)
m_k = m_k / 1000
D_k = ufloat(50.760, 0.004)
D_k = D_k / 1000
R_k = D_k / 2
Theta_halterung = 22.5
Theta_halterung = Theta_halterung / 10**7
L_draht = 0.585
Theta_kugel = 2 / 5 * m_k * R_k**2
print("Theta Kugel=", Theta_kugel)
Theta_ges = Theta_kugel + Theta_halterung
n, Ddraht = np.genfromtxt("Messdaten/daten_des_aufbau.txt", unpack=True)
Rdraht = Ddraht / 2
R_draht = np.mean(Rdraht)
R_draht = R_draht / 1000
m, t = np.genfromtxt("Messdaten/a.txt", unpack=True)
T = np.mean(t)

G = ((2 * L_draht) / (np.pi * R_draht**4)) * (2 * np.pi / T)**2 * Theta_ges
print("Schubmodul G=", G)
