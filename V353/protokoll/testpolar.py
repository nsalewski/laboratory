import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
#from astropy.io import ascii


w, urc, ug, a = np.genfromtxt("Messdaten/b_c.txt", unpack=True)

s = a / 1000
b = 1 / w
phi = 2 * np.pi * s / b
unull = 7.28
a_unull = urc / unull

RC = 3.67 * 10**(-3)
theta = np.linspace(0, 0.5 * np.pi)
m = np.linspace(1, 10000)
print(m*RC)
Aw =-(np.sin(theta) / m * RC)
plt.polar(theta, Aw, 'b-', label='Regression')
plt.polar(phi, a_unull, 'rx', label="Messwerte")
plt.savefig("build/testpolar.pdf")
