import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

n, v = np.genfromtxt("Messdaten/adrianundclemens/clemensadrian_v_null.txt", unpack=True)

ascii.write([n, v], 'Messdaten/v_null_ac.tex', format='latex')
