import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import sympy

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

    if err_var == None:
        err_vars == f.free_symbols

    for v in err_vars:
        err = Symbol('latex_std' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{ ' + latex(v) + '}'

    return latex(sympy.sqrt(s), symbol_names = latex_names)


t, A, B, a = sympy.var('t A B a')

f = (a * A * t^(a-1)) / (1 + B * t^a)^2
print(f)
print(error(f))
