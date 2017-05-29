import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
Z=[29,29,29,30,32,35,37,38,40,41,79,79]
Z=np.asarray(Z)
E=[8.048,8.028,8.905,9.65,11.1031,13.4737,15.1997,16.1046,17.9976,18.9856,13.7336,11.9187]
E=np.asarray(E)
E=E*1000
enull=1.6021766208*10**(-19)
R=13.6
r=10973731.568508
c=299792458
d=201.4*10**(-12)
h=ufloat(6.626070040*10**(-34),0.000000081*10**(-34))
Rzwo=r*h.nominal_value*c
R=Rzwo
E=E*enull
a=1/137
theta=(360/(2*np.pi))*(np.arcsin(((c*h.nominal_value)/(2*d*E))))
sigma=Z-np.sqrt(E/R-(a**2*Z**4)/4)
print(sigma)
print(theta)
ascii.write([Z,E/enull,np.round(theta,2),np.round(sigma,2)],'Tabelle_literatur.tex',format='latex')
grenz=5*(2*np.pi/360)
lambdi=2*d*np.sin(grenz)
print("lambda min= ",lambdi)
power=(h*c)/lambdi
power=power/enull
print("Emax= ",power)
Etheo=35*1000
Etheo=Etheo*enull
ltheo=(h*c)/Etheo
print("lamnda theo= ",ltheo )
xalpha1=[44.8,45.2]
yalpha1=[5.0,28.0]
xalpha2=[45.5,46.0]
yalpha2=[15.0,5.0]
xbeta1=[40,40.4]
ybeta1=[8.0,18.0]
xbeta2=[41.2,41.5]
ybeta2=[14.0,7.0]
def ausgabe(x,y,peak):
    m=np.sqrt(((y[0]-y[1])/(x[0]-x[1]))**2)
    b=y[0]-m*x[0]
    yval=peak/2
    value=(yval-b)/m
    value=value/2
    vtheta=value*(2*np.pi/360)
    Ehalb=(h*c)/(2*d*np.sin(vtheta))
    Ehalb=Ehalb/enull
    print(Ehalb)
    #print(m)
    #print(b)
    print("Wert= ",value)
betatheta=20.4*(2*np.pi/360)
alphatheta=22.6*(2*np.pi/360)
Ealpha=(h.nominal_value*c)/(2*d*np.sin(alphatheta))
Ebeta=(h.nominal_value*c)/(2*d*np.sin(betatheta))
sigmaone=29-np.sqrt(Ebeta/R-(a**2*29**4)/4)
sigmazwo=29-np.sqrt((Ebeta-Ealpha)/R-(a**2*29**4)/4)

print(sigmaone)
print(sigmazwo)

ausgabe(xalpha1,yalpha1,28)
ausgabe(xalpha2,yalpha2,28)
ausgabe(xbeta1,ybeta1,19)
ausgabe(xbeta2,ybeta2,19)
###########################################################################################
# Absorptionsspektren
Ry = ufloat(13.605693009, 0.000000084)

#Zink
E_zi = h*c/(2*d*np.sin(20.2*2*np.pi/360)*enull)
print('E_ZI = ', E_zi)
sigma_zn = 30-unp.sqrt(E_zi/Ry - (30**4)/(4*137**2))
print('sigma_zn= ', sigma_zn)

#Germanium
E_ge = h*c/(2*d*np.sin(16.3*2*np.pi/360)*enull)
print('E_GE = ', E_ge)
sigma_ge = 32-unp.sqrt(E_ge/Ry - (32**4)/(4*137**2))
print('sigma_ge= ', sigma_ge)

#Brom
E_br = h*c/(2*d*np.sin(13.35*2*np.pi/360)*enull)
print('E_BR = ', E_br)
sigma_br = 35-unp.sqrt(E_br/Ry - (35**4)/(4*137**2))
print('sigma_br= ', sigma_br)

#Zirkonium
E_zr = h*c/(2*d*np.sin(10*2*np.pi/360)*enull)
print('E_ZR = ', E_zr)
sigma_zr = 40-unp.sqrt(E_zr/Ry - (40**4)/(4*137**2))
print('sigma_zr= ', sigma_zr)

#Gold
E_au_beta = h*c/(2*d*np.sin(13.0*2*np.pi/360)*enull)
E_au_gamma = h*c/(2*d*np.sin(15.2*2*np.pi/360)*enull)
print('Gold', E_au_beta/enull, 'zweite:', E_au_gamma)
sigma_au = 79-unp.sqrt(4*137*unp.sqrt((E_au_beta-E_au_gamma)/Ry)-5*(E_au_beta-E_au_gamma)/Ry) * (1+19/(32*137**2)*(E_au_gamma-E_au_beta)/Ry)**(1/2)
print('SIGMAGOLD= ', sigma_au)

#########################################################################################
Ek=[8.91,10.97,13.33,17.73]
Ek=np.asarray(Ek)*1000
Z=[30,32,35,40]
Zls=np.linspace(29,41)
def theorie(x,m,b):
    return m*x+b
ascii.write([np.sqrt(Ek),Z],'tab_007.tex',format='latex',names=["wurzel e","Z"])
plt.plot(Z,np.sqrt(Ek), 'rx', label="Messwerte")
params, covariance = curve_fit(theorie,Z,np.sqrt(Ek))
errors = np.sqrt(np.diag(covariance))
ryd=ufloat(params[0],errors[0])
ryd=ryd**2
print('m= Rydbeck=',ryd)
print('b=',params[1],errors[1])
plt.plot(Zls, params[0]*Zls+params[1], 'b-', label="Lineare Regression")
plt.ylabel(r"$\sqrt{E_\mathrm{K}}$/$\sqrt{\si{\electronvolt}}$")
plt.xlabel(r"Kernladungszahl $Z$")
plt.xlim(29,41)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('R.pdf')
