import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
d, n = np.genfromtxt("Messdaten/spiegelverschiebung.txt", unpack=True)
uebersetzung=5.046
d=d/100
d=d/uebersetzung
wellenlaenge=(2*d)/n

ascii.write([np.round(d*10000,2),n, np.round(wellenlaenge * 10**9 , 2)], 'Messdaten/tab_wellenlaenge.tex', format="latex",
            names=['d','n','wellennlänge in nanometerz'])

wellen=ufloat(np.mean(wellenlaenge), np.std(wellenlaenge))
print("Wellenlänge: ", wellen)



p1, n1 = np.genfromtxt("Messdaten/luft.txt", unpack=True)
T0 = 273.15
p0 = 1.0132
T = 293.15
b = 50 * 10**(-3)
lambdi = 635 * 10**(-9)
ascii.write([p1, n1, np.round(1+T/T0 * p0/p1 * n1 * lambdi / (2*b), 6)], 'Messdaten/tab_lufti.tex', format = "latex")

luftibrechi = ufloat(np.mean(1+T/T0 * p0/p1 * n1 * lambdi / (2*b)), np.std(1+T/T0 * p0/p1 * n1 * lambdi / (2*b)))
print("Luftbrechi: ", luftibrechi)
print("Luftbrechi Δn: ", luftibrechi-1)


p2, n2 = np.genfromtxt("Messdaten/co2.txt", unpack=True)
ascii.write([p2, n2, np.round(1+T/T0 * p0/p2 * n2 * lambdi / (2*b), 6)], 'Messdaten/tab_c02.tex', format = "latex")

co2brechi = ufloat(np.mean(1+T/T0 * p0/p2 * n2 * lambdi / (2*b)), np.std(1+T/T0 * p0/p2 * n2 * lambdi / (2*b)))
print("co2brechi: ", co2brechi)
print("co2brechi Δn: ", co2brechi-1)
