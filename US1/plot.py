import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
pos1=1.374E+00
pos2=1.211E+02
print (pos2-pos1)
def form(t,c):
    return 0.5*c*t
def theorie(x,m,b):
    return m*x+b
s=120.5*10**(-3)
t=88.69*10**(-6)
c=2*s/t
print("Schallgeschwindigkeit=",c)
s,t=np.genfromtxt("Messdaten/laufzeit.txt", unpack=True)
ascii.write([s,t],'Messdaten/laufzeit.tex',format='latex',names=['s','t'])
s=s*10**(-3)
t=t*10**(-6)
params, covariance = curve_fit(theorie,t,s)
errors = np.sqrt(np.diag(covariance))
c_halbe=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
c=2*c_halbe
print("c= ",c)
print("b= ",b)
plt.plot(t*10**(6), s*10**3, 'rx', label="Messwerte")
t_l=np.linspace(20,90)
t_l=t_l*10**(-6)
plt.plot(t_l*10**(6), theorie(t_l,*params)*10**3, 'b-', label="Theoriekurve")

plt.ylabel(r"$s$/$\si{\milli\meter}")
plt.xlabel(r"$t$/$\si{\micro\second}$")
plt.xlim(20,90)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/ascan.pdf')
plt.clf()
s,t=np.genfromtxt("Messdaten/durchschall_messung_d.txt", unpack=True)
ascii.write([s,t],'Messdaten/durchschall.tex',format='latex',names=['s','t'])
s=s*10**(-3)
t=t*10**(-6)
params, covariance = curve_fit(theorie,t,s)
errors = np.sqrt(np.diag(covariance))
c=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
print("c Durchschall= ",c)
print("b Durchschall= ",b)
plt.plot(t*10**(6), s*10**3, 'rx', label="Messwerte")
t_ls=np.linspace(11,40)

t_ls=t_ls*10**(-6)
plt.plot(t_ls*10**(6), theorie(t_ls,*params)*10**3, 'b-', label="Theoriekurve")
plt.xlim(11,40)
plt.ylabel(r"$s$/$\si{\milli\meter}")
plt.xlabel(r"$t$/$\si{\micro\second}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/d.pdf')
plt.clf()
s,a1,a2=np.genfromtxt("Messdaten/amplitude.txt", unpack=True)
ascii.write([s,a1,a2,np.log(a1/a2)],'Messdaten/amplitude.tex',format='latex',names=['s','a1','a2','ln(a1/a2)'])

s=s*10**(-3)

params, covariance = curve_fit(theorie,s,np.log(a1/a2))
errors = np.sqrt(np.diag(covariance))
a=ufloat(params[0],errors[0])
b=ufloat(params[1],errors[1])
c=2*c_halbe
print("a= ",-a)
print("b= ",b)
s_ls=np.linspace(25,88)
s_ls=s_ls*10**(-3)
plt.plot(s*10**3, np.log(a1/a2), 'rx', label="Messwerte")
plt.plot(s_ls*10**3, theorie(s_ls,*params), 'b-', label="Theoriekurve")
plt.xlim(25,88)
plt.ylabel(r"$\ln\left(\frac{I(s)}{I_0}\right)$")
plt.xlabel(r"$s$/$\si{\milli\meter}$")

plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/amplitude.pdf')
plt.clf()
###############################################################################################
c_ = [2738, 2770]
t = [30.41*10**(-6), 34.87*10**(-6), 42.32*10**(-6)]
print('c = ', np.mean(c_), '+/-', np.std(c_))
c__ = unp.uarray(np.mean(c_), np.std(c_))
print('s = ', 1/2 * c * t[0])
print('s = ', 1/2 * c * t[1])
print('s = ', 1/2 * c * t[2])

###############################################################################################
#Bestimmung der Abstände im Auge:
c_L = 2500
c_GK = 1410
n,t=np.genfromtxt('Messdaten/auge.txt', unpack=True)
t=t*10**(-6)

s_1 = 1/2 * c_GK * (t[1]-t[0])
s_2 = 1/2 * c_GK * (t[2]-t[1])
s_3 = 1/2 * c_L * (t[3]-t[2]) + s_2
s_4 = 1/2 * c_GK * (t[4]-t[3]) + s_3
print('s_1= ', s_1)
print('s_2= ', s_2)
print('s_3= ', s_3)
print('s_4= ', s_4)
###############################################################################################

#n,f=np.genfromtxt("Messdaten/b_2.txt",unpack=True)
#f=f*1000
#theta=(n*np.pi)/14
#w=f*2*np.pi
#L=1.217*1/10**3
#C=20.13*1/10**9
#thetaplot = np.linspace(0, 3)
#def theorie(theta):
#    return np.sqrt(2/(L*C)*(1-np.cos(theta)))
#ascii.write([n,f/1000,np.round(f*2/1000*np.pi,1),np.round(theta,2)], 'Messdaten/tab_b1.tex', format="latex",
#            names=['n','frequenz','kreis','theta'])
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
