import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import sympy
from uncertainties import ufloat
import uncertainties.unumpy as unp
from sympy import Symbol, latex
from sympy import *

#aufgabe a: Messdaten plotten
t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60
T1=T1+273.15
T2=T2+273.15
pb=pb+1
pa=pa+1
plt.plot(t, T1, 'o', label='T1')
plt.plot(t, T2, 'o', label= 'T2')

#Aufgabe b: Curvefit
#f1+2
#def f1(t,a,b,c):
#    return a*t**2+b*t+c
#def f2(t,a,alpha,b):
#     return (a/(1+b*t**(alpha)))

#params1, covariance1 = curve_fit(f, t, T1)
#params2, covariance2 = curve_fit(f, t, T2)


#f3

def f3(t,a,alpha,b,c):
        return(((a*t**alpha)/(1+b*t**alpha)+c))
params1, covariance1 = curve_fit(f3, t, T1,bounds=([-np.inf,1,-np.inf, -np.inf],[np.inf,2,np.inf, np.inf]))
params2, covariance2 = curve_fit(f3, t, T2,bounds=([-np.inf,1,-np.inf, -np.inf],[np.inf,2,np.inf, np.inf]))

errors1 = np.sqrt(np.diag(covariance1))
errors2 = np.sqrt(np.diag(covariance2))

#hilfsfunktionen

ascii.write([t, T1, T2, pb, pa, leistung], 'hilfsdateien/values.dat', format='latex')
ascii.write([params1,errors1], 'hilfsdateien/paramserr1.dat', format='latex')
ascii.write([params2,errors2], 'hilfsdateien/paramserr2.dat', format='latex')


plt.plot(t,f3(t,*params1), label ="Fit für T1")
plt.plot(t,f3(t,*params2), label ="Fit für T2")
plt.grid()
plt.legend(loc='best')
plt.xlabel(r'$\frac{t}{s}$')
plt.ylabel(r'$Temperaturen \cdot\frac{1}{K}$')
plt.xlim(0, t[len(t)-1]+80)
plt.tight_layout()
plt.savefig("build/plot.pdf")


def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()

    if err_vars == None:
        err_vars = f.free_symbols

    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'

    return latex(sqrt(s), symbol_names=latex_names)

t, A, B, a = sympy.var('t A B a')

f = (a * A * t**(a-1)) / (1 + B * t**a)**2
print(f)
print(error(f))
#aufgabe e

A, B, a, L = sympy.var('A B a L')
c=17477
g = ((a * A * t**(a-1)) / (1 + B * t**a)**2)*1/L*c
A= ufloat(0.006,0.001)
a= ufloat(1.19,0.04)
B= ufloat(0.00011,0.00002)
C= ufloat(293.2,0.2)
L= ufloat(16770.3,124.7)
print("Massendurchsatz")
kappa=1.14
m=120.91

print(g)
print(error(g))
zeit=[300,600,900,1200]
def delta_m(a,A,B,t, L, c, m):
    return (((a * A * t**(a-1)) / (1 + B * t**a)**2)*(1/L)*c*m)
def rho(pa, T2):
    return (273.15*5.51*pa[t/60])/(T2[t/60])
for t in zeit:
    print(delta_m(a,A,B,t, L, c, m)/m)
print("Massendurchsazt in SI einheit")
for t in zeit:
    print(delta_m(a,A,B,t, L, c, m))
#print(g)
#print(error(g))


fehler=0
#aufgabe f
for t in zeit:
    fehler+=((1/(kappa-1))*(((pb[t/60])*((pa[t/60])/pb[t/60])**(1/kappa))-pa[t/60])*10**(2)*(1/rho(pa, T2))*(delta_m(a,A,B,t, L, c, m)))
    print ((1/(kappa-1))*(((pb[t/60])*((pa[t/60])/pb[t/60])**(1/kappa))-pa[t/60])*10**(2)*(1/rho(pa, T2))*(delta_m(a,A,B,t, L, c, m)))
    #*1000 für einheitenumrechnung watt=kg*m^2/s^3; rechte seite=bar*m^3/kg*g/s=m^2/s^2*g/s=>m^2/s^2*g*1000/s
print(np.mean(fehler))
for t in zeit:
    print(rho(pa, T2))
