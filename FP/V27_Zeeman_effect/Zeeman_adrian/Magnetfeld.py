import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from uncertainties import ufloat
from scipy.optimize import curve_fit

I, B =  np.genfromtxt("Magnetfeld.txt",unpack=True)
def f(a,b,c,d,x):
    return a*x**3+b*x**2+c*x+d


params,covariance= curve_fit(f,I,B)
errors=np.sqrt(np.diag(covariance))
x=np.linspace(0,20)

print("Magnetfeld bei 5 Ampere: ",f(5,*params))

plt.plot(I,B,"rx",label="Messwerte")
plt.plot(x,f(x,*params),'b--',label="Fit")
#print('Polynomparameter: a,b,c',*params,*errors)
plt.xlabel(r'Strom $ I/\si{\milli\A}$')
plt.ylabel(r'Magnetfeld $B/\si{\milli\tesla}$')
plt.tight_layout()
plt.legend(loc='best')
#plt.show()
plt.savefig("Eichung.pdf")
plt.clf()
