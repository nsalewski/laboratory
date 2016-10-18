import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


x=np.genfromtxt("1durchT.txt", unpack=True)
y=np.genfromtxt("ln(p).txt",unpack=True)
a=(np.mean(x))
b=(np.mean(y))
c=(np.mean(x*y))
d=(np.mean(x**2))
e=((np.mean(x))**2)
print(a,b,c,d,e)
print ((c-(a*b))/(d-e))
print(((d*b)-(a*c))/(d-e))
