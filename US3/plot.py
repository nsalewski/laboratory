import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

######
def theorie(x,m,b):
    return m*x+b
#########################################################
P ,f15, f30, f60, m15, m30, m60, d15, d30, d60=np.genfromtxt("Messdaten/a.txt", unpack=True)
# Delta v = 2 v0 v/c cos(a)
f15=np.absolute(f15)
f60=np.abs(f60)
m30=np.abs(m30)
d30=np.abs(d30)
v0 = 2*10**6
c = 1800
#Berechnet Dopplerwinkel
def local_doppler(theta):
    a=90-np.arcsin(np.sin(theta*2*np.pi/360)*2/3)*360/(2*np.pi)
    return a*2*np.pi/360
#print(local_doppler(15))

def local_pace(Delta, a):
    return (Delta * c/(2*v0*np.cos(a)))

#print(local_pace(f15, local_doppler(15)))
v0=2*10**(6)
c=1800
print('Steigung Theorie',2*v0/c)

#habs nur eben kurz programmiert, die labels passen alle noch nicht.
def local_pace_plot(fname,pname,f,m,b):#Funktion soll für jeden Prismenwinkel die nötigen Plots erstellen
    #Ausgleichgrade füür alle Rohrdicken
    deltanu=np.concatenate((f,m,b),axis=0)
    print(deltanu)
    params, covariance = curve_fit(theorie,local_pace((deltanu),local_doppler(pname)),(deltanu)/np.cos(local_doppler(pname)))
    errors = np.sqrt(np.diag(covariance))
    steigung=ufloat(params[0],errors[0])
    y=ufloat(params[1],errors[1])
    print('Steigung bei {} °'.format(pname), steigung)
    print('Achsenabschnitt bei {} °'.format(pname), y)

    plt.plot(local_pace((deltanu),local_doppler(pname)),(deltanu)/np.cos(local_doppler(pname)) , 'rx', label=r"Messwerte Röhre \theta={}\textdegree".format(pname))
    plt.plot(local_pace((deltanu),local_doppler(pname)),theorie(local_pace((deltanu),local_doppler(pname)),*params),'r-',label=r"Ausgleichgrade bei \theta={}\textdegree".format(pname))
    plt.ylabel(r"$\frac{\Delta \nu}{\cos{\alpha}}$/$\si{\Hz}$")
    plt.xlabel(r"$v$/$\si{\meter\per\second}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/Theta_{}.pdf'.format(fname))
    plt.clf()
local_pace_plot('15',15,f15,m15,d15)
local_pace_plot('60',60,f60,m60,d60)

local_pace_plot('30',30,f30,m30,d30)

#ploty = []
#for i in range(0, 4):
#    ploty[i]=f15[i]/np.cos(15)
#    ploty[i+3]=f30[i]/np.cos(30)
#    ploty[i+6]=f60[i]/np.cos(60)
#
#print(ploty)

########################################################
#Strömungsprofil

def local_plot(fname,pname):
    t,I,f=np.genfromtxt("Messdaten/b_{}.txt".format(fname),unpack=True)

    a=local_doppler(15)
    api=a*2*np.pi/360
    f0=2*10**6
    c=1800
    v=f/(2*f0)*c/np.cos(api)
    e=t-12.29
    e=e*1.798
    e=e+12.29*2.7
    ascii.write([t,np.round(e,1),I,f,np.round(v,3)],'Messdaten/tab_{}'.format(fname),format='latex',names=[r"Eindringtiefe in sec",r"Eindringtiefe/$\si{\milli\meter}$",r"$I_\mathrm{S}$/$\si{\square\volt\per\second}$",r"$\Delta \nu$/$\si{\Hz}$",r"momentane Geschwindigkeit/$\si{\meter\per\second}$"])
    t_ls=np.linspace(np.min(e)-0.1*(np.max(e)-np.min(e)),np.max(e)+0.1*(np.max(e)-np.min(e)))
    plt.plot(e, I, 'rx', label="Messwerte für P={}\%".format(pname))
    plt.ylabel(r"$I_\mathrm{S}$/$\si{\square\volt\per\second}$")
    plt.xlabel(r"$x/\si{\milli\meter}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/I{}.pdf'.format(fname))
    plt.clf()

    plt.plot(e, v, 'rx', label="Messwerte für P={}\%".format(pname))
    plt.ylabel(r"$v_\mathrm{mom}$/$\si{\meter\per\second}$")
    plt.xlabel(r"$x/\si{\milli\meter}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/f{}.pdf'.format(fname))
    plt.clf()
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
