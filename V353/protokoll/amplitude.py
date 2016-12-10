import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii



w, urc,ug,phi=np.genfromtxt("Messdaten/b_c.txt", unpack=True)
unull=728
urc=urc*100
def f(x, a):
    return (unull/(np.sqrt(1+(x**2)*(a**2))))
params, covariance = curve_fit(f,w,urc)
errors = np.sqrt(np.diag(covariance))
print('a = (RC) =', np.sqrt(params[0]), '±', np.sqrt(errors[0]))
print(params)
ascii.write([urc,w ], 'Messdaten/b.tex', format="latex")

temp=(f(w,*params))#used temp, bc without temp there was really freaky and wrong behaviour in matplotlib
plt.plot(w, urc, 'rx', label="Messwerte")
plt.plot(w, temp, 'b-', label='Ausgleichsgerade')
plt.xlim(4.24,10000)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$A(\omega)$/ $\si{\milli \volt}$")
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/amplitude.pdf")
