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
ascii.write([u,impsroh,np.round(imps,2),np.round((np.sqrt(imps)/imps),3)],'Messdaten/tab_a.tex',format='latex',names=["U in V","Imps mess","Imps per second","err"])
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

plt.errorbar(u,imps, yerr=[-np.sqrt(imps)/imps,+np.sqrt(imps)/imps], xerr=None, fmt='rx', label="Messwerte")
plt.ylabel(r"$N$/$\frac{\mathrm{[Imps]}}{\si{\second}}$")
plt.xlabel(r"$U$/$\si{\volt}$")
plt.xlim(min(u)-0.1*min(u),max(u)+0.1*max(u))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/a.pdf')
plt.clf()
print('*****************ende a************************')
u,t,te=np.genfromtxt('Messdaten/c.txt', unpack=True)
tm=np.mean(t)
ascii.write([u,t,te],'Messdaten/tab_c.tex', format='latex', names=['U','T','Te'])
#Wird nicht verwendet, da kein sinnvoller zusammenhang erkennbar
params, covariance = curve_fit(theorie,u, te)
errors = np.sqrt(np.diag(covariance))
m=ufloat(params[0],errors[0])
b=ufloat(params[1], errors[1])
print('Messung Erholungszeit: ')
print("m= ", m)
print("b= ", b)
uls=np.linspace(min(u)-0.1*min(u),max(u)+0.1*max(u))
plt.plot(uls,theorie(uls,*params), 'b-', label="Ausgleichsgrade")
plt.plot(u,te,'rx',label='Messwerte')
plt.ylabel(r"$T_\mathrm{E}$/$\si{\micro\second}$")
plt.xlabel(r"$U$/$\si{\volt}$")
plt.xlim(min(u)-0.1*min(u),max(u)+0.1*max(u))
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/c.pdf')
plt.clf()
print('***********************Ende Erholungszeit************************')
names=[r'N_1',r'N_{1+2}',r'N_{2}']
imps=np.asarray([3499,49676,46298])
imps_per_s=imps/60
imps1=ufloat(imps_per_s[0],np.sqrt(imps_per_s[0])/imps_per_s[0])
imps2=ufloat(imps_per_s[1],np.sqrt(imps_per_s[1])/imps_per_s[1])
imps3=ufloat(imps_per_s[2],np.sqrt(imps_per_s[2])/imps_per_s[2])
impstab=[imps1,imps2,imps3]

ascii.write([names,imps,impstab], 'Messdaten/totzeit.tex', format='latex')
T=(imps1+imps3-imps2)/(2*imps1*imps3)
print('Totzeit=',T)

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
