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
w,delta_t, uc, ut = np.genfromtxt("Messdaten/c_und_d.txt", unpack=True)

ascii.write([w,uc, np.round(uc/ut,2)], 'Messdaten/tab_c.dat', format="latex",
            names=[r"$f$/$\cdot 10^{3} \si{\Hz}$",r"$U_\text{C}$/$\si{\volt}$", r"$\frac{U(t)}{U_\text{C}}$"])
L=ufloat(10.11,0.03)
C=ufloat(2.098,0.006)
L=L/1000
C=C/10**9
R=ufloat(559.5,0.5)

f=w*1000
omega=R/L
omega=omega/(2*np.pi) #gotcha! that was the missing and huge errors causing line!
omeganull=unp.sqrt(1/(L*C))
q=1/(omeganull*R*C)
print("q=",q)
print("w-;w+=",omega)
plt.plot(f, uc/ut, 'rx', label="Messwerte")
plt.xlabel(r'$\omega$ /$\si{\Hz}$')
plt.ylabel(r"$\frac{U(t)}{U_\text{C}}$")
plt.xscale("log")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("build/c.pdf")
plt.clf()
x=[20,45]
y=[3.61/(np.sqrt(2)),3.61/(np.sqrt(2))]
x1=[28.5,28.5]
x2=[37.75,37.75]
y1=[1.5,4.0]

plt.plot(w[8:18], uc[8:18]/ut[8:18], 'rx', label="Messwerte")
plt.plot(x,y,"b--", label=r"$\frac{q_\mathrm{Experiment}}{\sqrt{2}}$")
plt.plot(x1,y1,"g--", label=r"$ \omega_\mathrm{+}$ und $\omega_\mathrm{-}$")
plt.plot(x2,y1,"g--")
plt.xlabel(r'$\omega$ /$\si{\Hz}$')
plt.ylabel(r"$\frac{U(t)}{U_\text{C}}$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("build/breite.pdf")
