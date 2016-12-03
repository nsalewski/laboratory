import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import sympy
from uncertainties import ufloat
import uncertainties.unumpy as unp
from sympy import Symbol, latex
from sympy import *
from pylab import *
ra=ufloat(15.5, 0.3)
ri=ufloat(10000000, 0)
u=ufloat(1.5,0.0225)
f=ra/ri
p=f/u
print("f=", f)
print("prozent=", p)
