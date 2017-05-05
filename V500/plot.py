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
    plt.ylim(np.min(np.sqrt(Itab)*10**5)-0.5,np.max(np.sqrt(Itab)*10**5)+1)
    plt.plot(x, lin(x,*params)*10**5, 'b-', label="linearer Fit")
    plt.errorbar(U, np.sqrt(Itab)*10**5, xerr=0, yerr=(errorI), fmt='rx', label="Messdaten samt Errorbalken")
    plt.ylabel(r"$\sqrt{I}\cdot 10^{-5}/\sqrt{\si{\ampere}}$")
    plt.xlabel(r"$U_{\mathrm{G}}/\si{\volt}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/{}.pdf'.format(filename))
    plt.clf()
    a=ufloat(params[0],errors[0])
    print('a von {}'.format(filename),a)
    b=ufloat(params[1],errors[1])
    print('b von {}'.format(filename),b)

    print('U von {}'.format(filename),-b/a)
    return -b/a

farben=["violett","blaugrün","grün","violett_drittkleinste","ultraviolett","gelb"]
Ug=[]
for x in farben:
    Ug.append(plot_and_return_ug(x))
U_G=np.asarray(Ug)
