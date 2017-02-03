import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
#

v = np.genfromtxt("Messdaten/a.txt", unpack="True")
v = v / 5
g = ufloat(9.811899, 0.000041)
l = ufloat(45, 0.3)
l = l / 100
links = v[0:10]
rechts = v[10:20]
gleichsinnig = v[20:30]
gegensinnig = v[30:40]
schwingung, schwebung = np.genfromtxt(
    "Messdaten/schwebunga.txt", unpack="True")
schwebung = schwebung / 5
ascii.write([links, rechts], 'Messdaten/einzelna.tex', format='latex')
ascii.write([gleichsinnig, gegensinnig],
            'Messdaten/koppelna.tex', format='latex')
ascii.write([schwingung, schwebung], 'Messdaten/schwebea.tex', format='latex')

links = ufloat(np.mean(links), np.std(
    links, ddof=1) / np.sqrt(len(links)))
rechts = ufloat(np.mean(rechts), np.std(
    rechts, ddof=1) / np.sqrt(len(rechts)))
gleichsinnig = ufloat(np.mean(gleichsinnig), np.std(
    gleichsinnig, ddof=1) / np.sqrt(len(gleichsinnig)))
gegensinnig = ufloat(np.mean(gegensinnig), np.std(
    gegensinnig, ddof=1) / np.sqrt(len(gegensinnig)))
schwebung = ufloat(np.mean(schwebung), np.std(
    schwebung, ddof=1) / np.sqrt(len(schwebung)))
schwingung = ufloat(np.mean(schwingung), np.std(
    schwingung, ddof=1) / np.sqrt(len(schwingung)))

print("Rechts", rechts)
print("Links", links)
print("Gleichsinnig", gleichsinnig)
print("Gegensinnig", gegensinnig)

omegaplus = ((2 * np.pi) / (gleichsinnig))
print("Omegaplus experiment", omegaplus)
omegaplustheo = unp.sqrt(g / l)
print("Omegaplus Theorie", omegaplustheo)
print("Tplus experimen", gleichsinnig)
tplustheo = (2 * np.pi * unp.sqrt(l / g))
print("Tplus Theorie", tplustheo)
kappa = (gleichsinnig**2 - gegensinnig**2) / (gleichsinnig**2 + gegensinnig**2)
print("Kopplungsgrad Kappa", kappa)
omegaminus = ((2 * np.pi) / (gegensinnig))
omegaminustheo = unp.sqrt((omegaplus**2 * (1 + kappa)) / (1 - kappa))
print("Omegaminus experiment", omegaminus)
print("Omegaminus theo", omegaminustheo)
print("Tminus experiment", gegensinnig)
tminustheo = 2 * np.pi * 1 / omegaminustheo
print("Tminus theorie", tminustheo)
omegas = omegaplus - omegaminus
print("Omega s Experiment", omegas)
omegastheorie = omegaplustheo - omegaminustheo
print("Omega s Theo", omegastheorie)


schwebungtheo = (tplustheo * tminustheo) / (tplustheo - tminustheo)
print("Schwebungsdauer theorie", schwebungtheo)
print("Schwingung", schwingung)

print("Schwebung", schwebung)
omegaex = 2 * np.pi / schwebung
omegatheo = 2 * np.pi / schwebungtheo
print("Omega s experiment über Zeit", omegaex)
print("Omega s theorie über Zeit", omegatheo)
print("****************************************************************************************")
