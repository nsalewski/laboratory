import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

T, p = np.genfromtxt('data2.txt', unpack = True)

def y(x, a, b, c,d):
    return a * x**3 + b*x**2 + c*x + d

T = T + 273.15


popt, cov = curve_fit(y, T, p)


plt.plot(T, p, 'rx',label='Messung')
l = np.linspace(330,451,1000)
plt.plot(l, y(l, *popt), 'b-', label='Fit')
plt.xlim(330, 455)
plt.xlabel(r'$T[K]$')
plt.ylabel(r'$p[bar] $')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messung2.pdf')

p=p*10000
popt, cov = curve_fit(y, T, p)
errors = np.sqrt(np.diag(cov))

print('Für die Koeffizienten gilt:')
print('a= ', popt[0], '+-', errors[0])
print('b= ', popt[1], '+-', errors[1])
print('c= ', popt[2], '+-', errors[2])
print('d= ', popt[3], '+-', errors[3])

R=8.3144621

def L(T,R,p,a,b,c,d,k):
    return  ((((R*T)/2) +p*  np.sqrt(     ( (R*T)/2)**2   -     (k* (a*T**3+b*T**2+c*T+d)) ))    *T* (3*a*T**2+b*2*T+c)/(a*T**3+b*T**2+c*T+d)  )

w=np.array([60, 80, 100, 120, 140, 160, 180])
w=w+273
o=np.array([42.48, 41.585, 40.657, 39.684, 38.643, 37.518, 36.304])
o=o*1000

plt.clf()

l=np.linspace(330,451)
plt.plot(l,L(l,R,1,popt[0],popt[1],popt[2],popt[3],0.9),'b-',label='$L_+$')
plt.plot(w,o,'rx',label='Literaturwerte')
plt.ylabel(r'$L(T)$')
plt.xlabel(r'$T[K] $')
plt.legend(loc='best')
plt.legend(loc='best')
plt.tight_layout()

plt.savefig('L+.pdf')

plt.clf()

plt.plot(l,L(l,R,-1,popt[0],popt[1],popt[2],popt[3],0.9),'b-',label='$L_-$')
plt.ylabel(r'$L(T)$')
plt.plot(w,o,'rx',label='Literaturwerte')
plt.xlabel(r'$T[K] $')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('L-.pdf')
