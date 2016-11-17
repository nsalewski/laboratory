import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import sympy
from uncertainties import ufloat
import uncertainties.unumpy as unp
from sympy import Symbol, latex
from sympy import *

r2,r3=np.genfromtxt("Messdaten/a.txt", unpack=True)
r4=1000-r3
rx=unp.uarray(r2*(r3/r4), r2*0.005*(r3*0.005/r4*0.005))
print(rx)
print(rx[0:3])
r_m=np.mean(rx[0:3])
print(r_m)
