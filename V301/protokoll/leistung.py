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

def y(x,Uo,ri):
	return (Uo**2*x)/(ri+x)**2

i,u = np.genfromtxt("Messdaten/b.txt", unpack=True)
u = u/1000
i = i/1000
ra = u/i
N = i**2 * ra
ri = 15.5
Uo = 1.5
x = np.linspace(0, 60, 500)

summi = 0
for i in range(0,len(ra)):
		summi = summi + N[i] - y(ra[i],Uo,ri)

n = len(ra)
summi = summi/len(ra)
sp = 0
for i in range(0,len(ra)):
		sp=sp+(N[i]-summi)**2

sp=np.sqrt((1/(n**2-n))*sp)
print('Das Mittel der Differenz : ', summi)
print('Mittlere Fehler : ' ,sp)
print(N[6])
print(y(26.1,Uo,ri))
plt.xlabel("Belastungswiderstand $R_{\\text{a}}$/$\\si{\\ohm}$")
plt.ylabel("abgegebende Leistung $N(R_{\\text{a}}$)/$\si{\\watt}$")
plt.plot(ra, N, 'rx', label="Messwerte")
plt.plot(x, y(x,Uo,ri),'b-',label="Theoriekurve")
plt.legend(loc="best")
plt.savefig("Bilder/leistung.pdf")
