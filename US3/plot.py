import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp







########################################################
#Strömungsprofil

def local_plot(fname,pname):
    t,I,f=np.genfromtxt("Messdaten/b_{}.txt".format(fname),unpack=True)
    t=t*3/2
    t_ls=np.linspace(np.min(t)-0.1*(np.max(t)-np.min(t)),np.max(t)+0.1*(np.max(t)-np.min(t)))
    plt.plot(t, I, 'rx', label="Messwerte für P={}\%".format(pname))
    plt.ylabel(r"$I_\mathrm{S}$/$\si{\square\volt\per\second}$")
    plt.xlabel(r"$x/\si{\milli\meter}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/I{}.pdf'.format(fname))
    plt.clf()

    plt.plot(t, f, 'rx', label="Messwerte für P={}\%".format(pname))
    plt.ylabel(r"$\Delta \nu$/$\si{\Hz}$")
    plt.xlabel(r"$x/\si{\milli\meter}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/f{}.pdf'.format(fname))
    plt.clf()
    ascii.write([t,I,f],'Messdaten/tab_{}'.format(fname),format='latex',names=[r"Eindringtiefe/$\si{\milli\meter}$",r"$I_\mathrm{S}$/$\si{\square\volt\per\second}$",r"$\Delta \nu$/$\si{\Hz}$"])
local_plot(70,70)
local_plot('45-2','45.2')

########################################################


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
