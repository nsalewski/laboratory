import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

def theorie(x,m,b):
    return m*x+b

impsroh, u, i=np.genfromtxt('Messdaten/a.txt',unpack=True)
imps=impsroh/60
impsreg=imps[4:len(imps)-3]
print(impsreg)
ascii.write([u,impsroh,np.round(imps,2)],'Messdaten/tab_a.tex',format='latex',names=["U in V","Imps mess","Imps per second"])
params, covariance = curve_fit(theorie,u[4:len(imps)-3], impsreg)
errors = np.sqrt(np.diag(covariance))
m=ufloat(params[0],errors[0])
b=ufloat(params[1], errors[1])
print('Messung a: ')
print("m= ", m)
print("b= ", b)
uls=np.linspace(min(u)-0.1*min(u),max(u)+0.1*max(u))
plt.plot(uls,theorie(uls,*params), 'b-', label="Ausgleichsgrade")
plt.axvline(x=480, ls=':', color="k", label="Arbeitsbereich des Zählrohrs")
plt.axvline(x=max(u)-30, ls=':', color="k")

plt.errorbar(u,imps, yerr=[-np.sqrt(imps),+np.sqrt(imps)], xerr=None, fmt='rx', label="Messwerte")
plt.ylabel(r"$N$/$\frac{\mathrm{[Imps]}}{\si{\second}}$")
plt.xlabel(r"$U$/$\si{\volt}$")
plt.xlim(min(u)-0.1*min(u),max(u)+0.1*max(u))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/a.pdf')
print('*****************************************')
#n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
#f=f*1000
#theta=(n*np.pi)/14
#w=f*2*np.pi
#L=1.217*1/10**3
#C=20.13*1/10**9
#thetaplot = np.linspace(0, 3)
#
#def theorie(theta):
#    return np.sqrt(2/(L*C)*(1-np.cos(theta)))
#
#ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
#            names=['n','frequenz','kreis','theta'])
#
#
#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#curve_fitting:
#def theorie(x,m,b):
#    return m*x+b
#ascii.write([np.sqrt(Ek),Z],'tab_007.tex',format='latex',names=["wurzel e","Z"])
#plt.plot(Z,np.sqrt(Ek), 'rx', label="Messwerte")
#params, covariance = curve_fit(theorie,Z,np.sqrt(Ek))
#errors = np.sqrt(np.diag(covariance))
#ryd=ufloat(params[0],errors[0])
