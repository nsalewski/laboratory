import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
d, n = np.genfromtxt("Messdaten/spiegelverschiebung.txt", unpack=True)
uebersetzung=5.024
d=d/100
d=d/uebersetzung
wellenlaenge=(2*d)/n

ascii.write([np.round(d*10000,2),n, np.round(wellenlaenge * 10**9 , 2)], 'Messdaten/tab_wellenlaenge.tex', format="latex",
            names=['d','n','wellennlänge in nanometerz'])

wellen=ufloat(np.mean(wellenlaenge), np.std(wellenlaenge))
print("Wellenlänge: ", wellen)
