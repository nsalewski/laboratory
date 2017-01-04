import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

n,f=np.genfromtxt("Messdaten/a_2.txt",unpack=True)
n=n*4

ascii.write([n,f], 'Messdaten/tab_a1.tex', format="latex",
            names=['n','frequenz'])


