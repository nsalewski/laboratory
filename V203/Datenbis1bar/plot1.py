import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

x1=np.linspace(0.0035,0.00265, 77)
x=np.genfromtxt("1durchT.txt", unpack=True)
y=np.genfromtxt("ln(p).txt",unpack=True)
a=(np.mean(x))
b=(np.mean(y))
c=(np.mean(x*y))
d=(np.mean(x**2))
e=((np.mean(x))**2)
print(a,b,c,d,e)
f=((c-(a*b))/(d-e))
g=((d*b)-(a*c))/(d-e)
print (f,g)
print (len(y))


plt.plot (x, y,"rx", label="Messdaten")
plt.plot(x1, f*x1+g, "b-", label="Regression")
plt.xlim(0.00265, 0.0034)

plt.grid()

plt.xlabel(r"$\frac{1}{T}\  in\  \frac{1}{K}$")
plt.ylabel("$log(p)$")
plt.title("Logarithmus des Drucks aufgetragen gegen den Kehrwert \n der Temperatur mit Regressionsgeraden")

plt.savefig("messung1.pdf")
