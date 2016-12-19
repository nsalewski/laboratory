import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
w, a, uc, u = np.genfromtxt("Messdaten/c_und_d.txt", unpack=True)
b = 1 / (w * 1000)
phi = 2 * np.pi * a / (10**6 * b)
#R = 559.5 #müsste doch 559.5 sein, wegen generatorinnenwiderstand. habs mal geändert. Vorher war es 509.5
L=ufloat(10.11,0.03)
C=ufloat(2.098,0.006)
L=L/1000
C=C/10**9
R=ufloat(559.5,0.5)
c=(R/(2*L))
print(c)
w1=(c +unp.sqrt((1/(L*C))+((R**2)/(4*L**2))))
w2=((-c) +unp.sqrt((1/(L*C))+((R**2)/(4*L**2))))
wres=unp.sqrt((1/(L*C))-((R**2)/(2*L**2)))
#quite not sure 'bout the following. Thought, it was already in correct frequency-type ()
w1=w1/(2*np.pi)
w2=w2/(2*np.pi)
wres=wres/(2*np.pi)


print("wres= ",wres)
print("w1=",w1)
print("w2=",w2)

def f(w):
    return np.arctan(-(w * R * C) / (1 - L * C * w**2))

ascii.write([w, a, b*10**6, phi], 'Messdaten/d.tex', format="latex")

#m = np.logspace(0.01, 4)
#temp = (f(m, *params))

plt.figure(0)
plt.plot(w, phi, 'rx', label="Messwerte")
#plt.plot(m, temp, 'b-', label='Ausgleichskurve')
plt.xlim(1, 100)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\phi(\omega)$/ $\si{\radian}$")
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/taskd.pdf")


#changed to manual lines, because its easier to control colours

x1=[10,45]
y1=[0.25*np.pi,0.25*np.pi]
y2=[0.5*np.pi,0.5*np.pi]
y3=[0.75*np.pi,0.75*np.pi]
y4=[0,5]
x2=[30.2,30.2]
x3=[33.7,33.7]
x4=[37.7,37.7]
plt.figure(1)
plt.plot(w, phi, 'rx', label="Messwerte")
plt.plot(x1,y1, "b-.")
plt.plot(x1,y2, "g-.")
plt.plot(x1,y3, "b-.")

plt.plot(x2,y4,"b--")
plt.plot(x3,y4,"g--",label="Resonatorfrequenz")
plt.plot(x4,y4,"b--",label="$\omega_1$ bzw.$\omega_2$")
plt.yticks( [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi],
            [ r'$0$', r'$\pi/4$',r"$\pi/2$" ,r'$3\pi/4$',r"\pi"]
    )
plt.xlim(10, 45)
plt.ylim(0,np.pi)
plt.xlabel("$\omega$ / $\si{\Hz}$")
plt.ylabel(r"$\phi(\omega)$/ $\si{\radian}$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("build/taskdlinear.pdf")
