import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp

nr,d=np.genfromtxt("Messdaten/abmessungen.txt",unpack=True)
ascii.write([nr,d*10],'Messdaten/tab_abmessung.tex',format='latex', names=['$n$', '$d$/'r"$\si{\milli\meter}$"])

#######################################################################################
c=2730
number, side1, side2 = np.genfromtxt('Messdaten/ascan.txt', unpack = True)
side1 = side1 - 2
side2 = side2 - 2
position1 = 1/2 * c * side1 * 10**(-6)
position2 = 1/2 * c * side2 * 10**(-6)
ascii.write([number, side1, position1, side2, position2],'Messdaten/tab_.ascantex',format='latex')
ascii.write([number, 0.08 - position1 - position2], 'Messdaten/tab_dicken.tex', format='latex')



##########Auswertung bscan###############
nr=np.linspace(1,11,11)
rueck=np.asarray([489,500,136,199,265,330,392,453,514,2, 149])
hin=np.asarray([-1,166,501,443,384,328,268,205,145,84,456])
rueck=rueck-22
hin=hin-22
pixel=724
sec=69
umrech=sec/pixel
def umrechnung(a,c):
    return c*a
hin=umrechnung(hin,umrech)
rueck=umrechnung(rueck,umrech)
anpass=0.9
hin=hin-anpass
rueck=rueck-anpass
t_rueck=rueck
t_hin=hin
hinweg=0.5*c*hin*10**(-4)
rueckweg=0.5*c*rueck*10**(-4)
ascii.write([nr,np.round(t_hin,1),np.round(t_rueck,1),np.round(hinweg,2),np.round(rueckweg,2),np.round(8-hinweg-rueckweg,2),d,np.round((np.round(8-hinweg-rueckweg,2)-d)/d*100,1)],'bscan/tab_bscan.tex',format='latex',names=['$n$','thin','trueck','hin', 'rueck','d','dnom','diff'])

################Auswertung Herz################
maximum=np.asarray([430,400,399,391,414,396,392,395,396,407,398,399,396])
minimum=np.asarray([535,540,536,538,534,535,535,533,535,538,537,535,539])
maximum=maximum-22
minimum=minimum-22
umrech2=69/706
minimum=umrechnung(minimum,umrech2)
maximum=umrechnung(maximum,umrech2)
c=1497
minweg=0.5*c*minimum*10**(-4)
maxweg=0.5*c*maximum*10**(-4)
ascii.write([np.round(minimum,2),np.round(maximum,2),np.round(minweg,2),np.round(maxweg,2)],'bscan/herz.tex',format='latex',names=['min','max','hmin','hmax'])
medmin=np.mean(minweg)
medmax=np.mean(maxweg)
print('hmin, hmax=',medmin-medmax)
r=6
h=medmin-medmax
V=np.pi/3*h**2*(3*r-h)
f=1.3
print('Schlagvolumen= ',V)
print('HZV=',V*f)

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
