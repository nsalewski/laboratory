import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


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
#

def lin(x,m,b):
    return m*x+b

Uvio,Ivio=np.genfromtxt("Messdaten/violett.txt",unpack=True)
Iviotab=Ivio/10**9
ascii.write([Uvio,Ivio, np.sqrt(Iviotab)*10**5],"Messdaten/tab_violett.tex", format="latex",names=["U","I*10^(-9)","sqrt(I)*10^(-5)"])
paramsvio, covariancevio = curve_fit(lin,Uvio ,np.sqrt(Iviotab))
errors1 = np.sqrt(np.diag(covariancevio))
errYvio = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
errYvio=[((np.sqrt(x/10**9))*10**5)for x in errYvio]
x=np.linspace(-0.2,2.2)
plt.ylim(-0.5,7.8 )
plt.xlim(-0.2,1.6)
plt.plot(x, lin(x,*paramsvio)*10**5, 'b-', label="linearer Fit")
plt.errorbar(Uvio, np.sqrt(Iviotab)*10**5, xerr=0, yerr=(errYvio), fmt='rx', label="Messdaten samt Errorbalken")
plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}/\sqrt{\si{\ampere}}$")
plt.xlabel(r"$U_{\mathrm{G}}/\si{\volt}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/violett.pdf')
plt.clf()

Uvioo,Ivioo=np.genfromtxt("Messdaten/ultraviolett.txt",unpack=True)
Iviootab=Ivioo/10**9
ascii.write([Uvioo,Ivioo, np.sqrt(Iviootab)*10**5],"Messdaten/tab_ultraviolett.tex", format="latex",names=["U","I*10^(-9)","sqrt(I)*10^(-5)"])
paramsvioo, covariancevioo = curve_fit(lin,Uvioo ,np.sqrt(Iviootab))
errors2 = np.sqrt(np.diag(covariancevioo))
errYvioo = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
errYvioo=[((np.sqrt(x/10**9))*10**5)for x in errYvioo]
x=np.linspace(-0.2,2.2)
plt.ylim(-0.5,10 )
plt.xlim(-0.2,2.0)
plt.plot(x, lin(x,*paramsvioo)*10**5, 'b-', label="linearer Fit")
plt.errorbar(Uvioo, np.sqrt(Iviootab)*10**5, xerr=0, yerr=(errYvioo), fmt='rx', label="Messdaten samt Errorbalken")
plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}/\sqrt{\si{\ampere}}$")
plt.xlabel(r"$U_{\mathrm{G}}/\si{\volt}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/ultraviolett.pdf')
plt.clf()


Uviu,Iviu=np.genfromtxt("Messdaten/violett_drittkleinste.txt",unpack=True)
Iviutab=Iviu/10**9
ascii.write([Uviu,Iviu, np.sqrt(Iviutab)],"Messdaten/tab_violett3.tex", format="latex",names=["U","I*10^(-9)","sqrt(I)*10^(-5)"])
paramsviu, covarianceviu = curve_fit(lin,Uviu ,np.sqrt(Iviutab)*10**5)
errors3 = np.sqrt(np.diag(covarianceviu))
errYviu = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
errYviu=[((np.sqrt(x/10**9))*10**5)for x in errYviu]
x=np.linspace(-0.2,2.2)
plt.ylim(-0.5,9.2)
plt.xlim(-0.1,1.4)
plt.plot(x, lin(x,*paramsviu)*10**5, 'b-', label="linearer Fit")
plt.errorbar(Uviu, np.sqrt(Iviutab)*10**5, xerr=0, yerr=(errYviu), fmt='rx', label="Messdaten samt Errorbalken")
plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}/\sqrt{\si{\ampere}}$")
plt.xlabel(r"$U_{\mathrm{G}}/\si{\volt}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/violett3.pdf')
plt.clf()

Uvib,Ivib=np.genfromtxt("Messdaten/blaugrün.txt",unpack=True)
Ivibtab=Ivib/10**9
ascii.write([Uvib,Ivib, np.sqrt(Ivibtab)*10**5],"Messdaten/tab_blaugrün.tex", format="latex",names=["U","I*10^(-9)","sqrt(I)*10^(-5)"])
paramsvib, covariancevib = curve_fit(lin,Uvib ,np.sqrt(Ivibtab))
errors4 = np.sqrt(np.diag(covariancevib))
errYvib = [0.01,0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
errYvib=[((np.sqrt(x/10**9))*10**5)for x in errYvib]
x=np.linspace(-0.2,2.2)
plt.ylim(-0.5,2.6)
plt.xlim(-0.1,1.4)
plt.plot(x, lin(x,*paramsvib)*10**5, 'b-', label="linearer Fit")
plt.errorbar(Uvib, np.sqrt(Ivibtab)*10**5, xerr=0, yerr=(errYvib), fmt='rx', label="Messdaten samt Errorbalken")
plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}/\sqrt{\si{\ampere}}$")
plt.xlabel(r"$U_{\mathrm{G}}/\si{\volt}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/blaugrün.pdf')
plt.clf()
