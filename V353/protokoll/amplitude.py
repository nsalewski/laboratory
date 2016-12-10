import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit



w, urc,ug,phi=np.genfromtxt("Messdaten/b_c.txt", unpack=True)
unull=7.28
def f(x, a):
    return (unull)/(np.sqrt(1+x**2*a**2))
params, covariance = curve_fit(f,w,urc)
errors = np.sqrt(np.diag(covariance))
print('a **2 = (RC)**2 =', params[0], '±', errors[0])
print(params)
m=np.linspace(0,10000)
plt.plot(w, urc/unull, 'rx', label="Messwerte")
plt.plot(m, f(m, *params), 'b-', label='Ausgleichsgerade')
#plt.xlabel(blaa)
#plt.ylabel(blup)
plt.legend(loc='best')
plt.savefig("build/amplitude.pdf")
