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
f15=np.abs(f15)
f60=np.abs(f60)
m30=np.abs(m30)
d30=np.abs(d30)
v0 = 2*10**6
c = 1800
#Berechnet Dopplerwinkel
def local_doppler(theta):
    a=90-np.arcsin(np.sin(theta*2*np.pi/360)*2/3)*360/(2*np.pi)
    return a*2*np.pi/360
v0=2*10**(6)
c=1800
print('Steigung Theorie',2*v0/c)

def vP(d, P):
    return P/100 * 0.01/((d/2)**2*np.pi)/60

v_f = np.array(vP(0.016, P))
v_m = np.array(vP(0.010, P))
v_d = np.array(vP(0.007, P))
print(v_f,v_m,v_d)
def local_plot(fname,v,nu,winkel): #Verwendung fname ist dateiname, v ist geschwindigkeit der Flüssigkeit
#nu ist momentane Dopplerverschiebung, winkel ist zugehöriger theta-Winkel
    params, covariance = curve_fit(theorie,v,(nu)/np.cos(local_doppler(winkel)))
    errors = np.sqrt(np.diag(covariance))
    steigung=ufloat(params[0],errors[0])
    y=ufloat(params[1],errors[1])
    print('Steigung bei {} {}°'.format(fname,winkel), steigung)
    print('Achsenabschnitt {} {}°'.format(fname,winkel), y)
    plt.plot(v,(nu)/np.cos(local_doppler(winkel)) , 'rx', label=r"Messwerte {} {}\textdegree".format(fname,winkel))
    plt.plot(v,theorie(v,*params),'b-',label=r"Ausgleichgrade {} {}\textdegree".format(fname,winkel))
    plt.ylabel(r"$\frac{\Delta \nu}{\cos{\alpha}}$/$\si{\Hz}$")
    plt.xlabel(r"$v$/$\si{\meter\per\second}$")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Bilder/{}_{}.pdf'.format(fname,winkel))
    plt.clf()
    return(steigung)#gibt steigung zum mittelwertbilden zurück

#Arrays zum einfacheren Funktionsaufruf
v_arr=np.asarray([v_f,v_f,v_f,v_m,v_m,v_m,v_d,v_d,v_d])
nu_arr=np.asarray([f15,f30,f60,m15,m30,m60,d15,d30,d60])
namearr=np.asarray(['dick','dick','dick','mittel','mittel','mittel','dünn','dünn','dünn'])
degree_arr=np.asarray([15,30,60,15,30,60,15,30,60])

steigung_arr=[]
#funktionenaufruf und steigung in ein array speichern
for i in range(0,9):
    steigung_arr.append(local_plot(namearr[i],v_arr[i],nu_arr[i],degree_arr[i]))

steigung=np.asarray(steigung_arr)
steigung=np.mean(steigung)
print(steigung)
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


#Reste von 1:
# Delta v = 2 v0 v/c cos(a)
#print(local_doppler(15))

#def local_pace(Delta, a):
#    return (Delta * c/(2*v0*np.cos(a)))

#print(local_pace(f15, local_doppler(15)))

#def local_pace_plot(fname,pname,v,f,m,b):#Funktion soll für jeden Prismenwinkel die nötigen Plots erstellen
#    print(v)
#    deltanu=np.concatenate((f,m,b),axis=0)
#    params, covariance = curve_fit(theorie,v,(deltanu)/np.cos(local_doppler(pname)))
#    errors = np.sqrt(np.diag(covariance))
#    steigung=ufloat(params[0],errors[0])
#    y=ufloat(params[1],errors[1])
#    print('Steigung bei {} °'.format(pname), steigung)
#    print('Achsenabschnitt bei {} °'.format(pname), y)
#    plt.plot(v,(deltanu)/np.cos(local_doppler(pname)) , 'rx', label=r"Messwerte Röhre \theta={}\textdegree".format(pname))
#    #plt.plot(local_pace((deltanu),local_doppler(pname)),(deltanu)/np.cos(local_doppler(pname)) , 'rx', label=r"Messwerte Röhre \theta={}\textdegree".format(pname))
#    plt.plot(v,theorie(v,*params),'r-',label=r"Ausgleichgrade bei \theta={}\textdegree".format(pname))
#    plt.ylabel(r"$\frac{\Delta \nu}{\cos{\alpha}}$/$\si{\Hz}$")
#    plt.xlabel(r"$v$/$\si{\meter\per\second}$")
#    plt.legend(loc='best')
#    plt.tight_layout()
#    plt.savefig('Bilder/Theta_{}.pdf'.format(fname))
#    plt.clf()
#local_pace_plot('15',15,vall,d15,m15,f15)
#local_pace_plot('30',30,vall,d30,m30,f30)
#local_pace_plot('60',60,vall,d60,m60,f60)
