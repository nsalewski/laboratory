import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#aufgabe a: Messdaten plotten
t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60

plt.plot(t, T1, 'o', label='Messdaten T1 (warmes Reservoir)')
plt.plot(t, T2, 'o', label= 'Messdaten T2 (kaltes Reservoir')

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
print(params1)
print(params2)



errors1 = np.sqrt(np.diag(covariance1))
errors2 = np.sqrt(np.diag(covariance2))

plt.plot(t,f3(t,*params1), label ="Fit für T1")
plt.plot(t,f3(t,*params2), label ="Fit für T2")

plt.savefig("build/plot.pdf")
