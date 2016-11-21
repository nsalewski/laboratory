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

r2,r3=np.genfromtxt("Messdaten/a.txt", unpack=True)
r4=1000-r3
rx=unp.uarray(r2*(r3/r4), r2*0.005*(r3*0.005/r4*0.005))
#print(rx)
#print(rx[0:3])
r_m=np.mean(rx[0:3])
#print(r_m)
f, us,ub=np.genfromtxt("Messdaten/e.txt", unpack=True)
c=658*10**-9
r=1000
w=2*np.pi*f*r*c
w0=1/r*c
v=1/9
ug=ub/2*np.sqrt(2)
ug=ug/1000
uc=ub/4000*np.sqrt(2)
um=ub/2*np.sqrt(2)
#print(ub)
g=logspace(np.log10(20), np.log10(30000), 31)
gamma=g/241.6

ubus=uc/us
ubus2=v*((gamma**2-1)**2/((1-gamma**2)**2+9*gamma**2))
print(ubus, ubus2)
ascii.write([f, us.round(2), um.round(2), (uc/us).round(4),np.sqrt(ubus2).round(4) ,(f/241.6).round(2)], "Messdaten/values.dat", format="latex")
#print(f)
plt.plot(f/241.6, ubus, 'o', label='Messdaten')
plt.xlabel(r'$\Omega= \frac{f}{f_0}$')
plt.ylabel(r"$\frac{U_{Br, eff}}{U_s}$")
plt.xscale("log")
plt.plot(g/241.6,np.sqrt(ubus2), '-', label= "Theoriekurve")
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.show()

plt.savefig("build/plot.pdf")
