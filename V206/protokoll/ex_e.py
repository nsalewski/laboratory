import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii

#aufgabe e
t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60
T1=T1+273.15
T2=T2+273.15
pb=pb+1
pa=pa+1
T1=np.delete(T1, 0)
pb=np.delete(pb,0)
x=1/T1
y=np.log(pb)

def func(x,a,b):
    return a*x+b
params_e, covariance_e=curve_fit(func, x, y)
plt.plot (x, y,"rx", label="Messdaten")
plt.plot(x, func(x,*params_e), "b-", label="Regression")
plt.grid()
errors_e=np.sqrt(np.diag(covariance_e))
ascii.write([params_e, errors_e], 'hilfsdateien/exercise.dat', format='latex')

plt.xlabel(r'$\frac{1}{T_1} \cdot K$')
plt.ylabel(r"$\log(\frac{p_b}{p_0})$")
plt.title("Verdampfungswärme")

plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/ex_e.pdf")
