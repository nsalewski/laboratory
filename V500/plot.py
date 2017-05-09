import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


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

def plot_and_return_ug(filename):
    U,I=np.genfromtxt('Messdaten/{}.txt'.format(filename),unpack=True)
    Itab=I/10**9
    helplist=[]
    for x in I:
        if x<=1:
            helplist.append(0.01)
        else: helplist.append(0.1)
    errorI=np.asarray(helplist)
    ErrorI=[(x/10**9) for x in errorI]
    errorI=[(np.sqrt(x))*10**5 for x in ErrorI]
    uncertainty_array_I= unp.uarray(Itab,ErrorI)
    ascii.write([U,uncertainty_array_I*10**9,unp.sqrt(uncertainty_array_I)*10**5],'Messdaten/{}.tex'.format(filename),format="latex",names=["U","$I*10^(-9)$","$sqrt(I)*10^(-5)$"])
    if(filename=='gelb'):
        U=U[29:]
        Itab=Itab[29:]
        errorI=errorI[29:]
    elif (filename=='grün'):
        U=U[3:]
        Itab=Itab[3:]
        errorI=errorI[3:]
    x=np.linspace(np.min(U)-(np.max(U)-np.min(U))*0.15,np.max(U)+(np.max(U)-np.min(U))*0.15)
    params, covariance=curve_fit(lin,U,np.sqrt(Itab))
    errors=np.sqrt(np.diag(covariance))
    plt.xlim(np.min(U)-(np.max(U)-np.min(U))*0.15,np.max(U)+(np.max(U)-np.min(U))*0.15)
    plt.ylim(np.min(np.sqrt(Itab)*10**5)-0.6,np.max(np.sqrt(Itab)*10**5)+1.2)
    plt.plot(x, lin(x,*params)*10**5, 'b-', label="lineare Regressionsgrade")
    plt.errorbar(U, np.sqrt(Itab)*10**5, xerr=0, yerr=(errorI), fmt='ro',capsize=3, label="Messdaten samt Errorbalken")
    plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}$/$\sqrt{\si{\ampere}}$")
    plt.xlabel(r"$U_{\mathrm{B}}$/$\si{\volt}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/{}.pdf'.format(filename))
    plt.clf()
    a=ufloat(params[0],errors[0])
    print('a von {}'.format(filename),a)
    b=ufloat(params[1],errors[1])
    print('b von {}'.format(filename),b)

    print('U von {}'.format(filename),-b/a)
    bdiva=-params[1]/params[0]
    return (bdiva)

farben=["ultraviolett","violett","violett_drittkleinste","blaugrün","grün","gelb"]
Ug=[]
for x in farben:
    Ug.append(plot_and_return_ug(x))
U_G=np.asarray(Ug)
Wellenlänge=[365,405,435,492,546,577]
Wellenlänge=[(x/10**9) for x in Wellenlänge]
c=299792458.0
Frequenz=[(c/x) for x in Wellenlänge]
paramsg, covarianceg=curve_fit(lin,Frequenz,U_G)
errorsg=np.sqrt(np.diag(covarianceg))
a=ufloat(paramsg[0],errorsg[0])
b=ufloat(paramsg[1],errorsg[1])
print("a=h/e: ",a)
print("b=ak/e: ",-b)
Frequenztab=[(x/10**14) for x in Frequenz]
ascii.write([Wellenlänge,Frequenztab,U_G],'Messdaten/Ug.tex',format="latex",names=["lamnda","nü","Ug"])
params, covariance=curve_fit(lin,Frequenztab,U_G)
errors=np.sqrt(np.diag(covariance))
x=np.linspace(np.min(Frequenztab)-(np.max(Frequenztab)-np.min(Frequenztab))*0.15,np.max(Frequenztab)+(np.max(Frequenztab)-np.min(Frequenztab))*0.15)
plt.xlim(np.min(Frequenztab)-(np.max(Frequenztab)-np.min(Frequenztab))*0.15,np.max(Frequenztab)+(np.max(Frequenztab)-np.min(Frequenztab))*0.15)
plt.plot(x, lin(x,*params), 'b-', label="lineare Regressionsgrade")
plt.plot(Frequenztab,U_G,'rx',label="Messdaten")
plt.ylabel(r"$U_\mathrm{G}$/$\si{\volt}$")
plt.xlabel(r"$\nu \cdot 10^{14}$/$\si{\hertz}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/Ug.pdf')
plt.clf()
h=ufloat(6.626070040,0.000000081)
h=h/10**34
e=ufloat(1.6021766208, 0.0000000098)
e=e/10**19
he=h/e
print("Theorie h/e= ",he)

#######################################################################################


U_g, I_ = np.genfromtxt("Messdaten/gelb.txt", unpack = True)

plt.plot(U_g, I_, 'rx', label="Messdaten")
#plt.plot(U_g,  U_g*2, 'b-', label="Lessdaten")
plt.xlabel(r"$U$ / $\si{\volt}$")
plt.ylabel(r"$I$ / $\si{\nano\ampere}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/gelbplot.pdf')
plt.clf()
