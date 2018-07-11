#!usr/bin/env python
#coding:utf8
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
import scipy.constants as const
import math as math
x=np.linspace(0,50,50)
y=np.linspace(0,50,50)
for i in range(len(x)-2):
    plt.plot([x[i], x[i+1]],[y[i],y[i]], "b-")
plt.savefig("test.pdf")
