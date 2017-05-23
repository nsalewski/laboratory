import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
#Berechnung mittlere freie Weglänge
#die ersten beiden sind deine Plots.
T=[415.95, 423.15, 299.75, 458.9]
T=np.asarray(T)
psaet=5.5*10**7*np.exp(-6876/T)
w=0.0029/psaet
a=1
#print(w)
ascii.write([T,psaet,w,w/a],'Messdaten/w.tex',format='latex' )

n=np.genfromtxt("Messdaten/datenraumtemperatur.txt",unpack=True)
#print(n)
c=np.mean(n)
m=np.std(n) / np.sqrt(len(n))
v1=ufloat(c,m)
v1=1/v1
print('Skala Raumtemperatur',v1)
z,n,pos=np.genfromtxt("Messdaten/steigung.txt",unpack=True)
steigung=z/n
pos=pos*v1.nominal_value
ascii.write([np.round(pos,1),np.round(steigung,2)],'Messdaten/steigii.tex',format='latex')
#print(steigung)
#print(pos)
plt.plot(pos, steigung, 'rx', label="differentielle Energieverteilung")
plt.ylabel(r"prop. zur Steigung")
plt.xlabel(r"$U_\mathrm{A}$/$\si{\volt}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Bilder/b1.pdf')
plt.clf()



p=np.genfromtxt("Messdaten/skalafranckhertz.txt",unpack=True)
c=np.mean(p)
m_2=np.std(p) / np.sqrt(len(p))
vhertz=ufloat(c,m_2)
print('Skala Franckhertz',vhertz)

vhertz=5/vhertz
print('Skala Franckhertz',vhertz)
e_theo = unp.uarray(1.6021766208*10**(-19), 0.0000000098*10**(-19))
abstand=np.genfromtxt("Messdaten/maximafranckhertz.txt", unpack=True)
ascii.write([abstand,abstand*vhertz.nominal_value],'Messdaten/franck.tex',format='latex')
c=np.mean(abstand*vhertz.nominal_value)
m_3=np.std(abstand*vhertz.nominal_value)/np.sqrt(len(abstand*vhertz.nominal_value))
print(c,m_3)
d=np.mean(abstand)
m_4=np.std(abstand)/np.sqrt(len(abstand))
abst=ufloat(d,m_4)
delta=ufloat(c,m_3)
e_theo = unp.uarray(1.6021766208*10**(-19), 0.0000000098*10**(-19))
max1=44-2*abst.nominal_value
max1*vhertz
delta=delta*e_theo
print(max1)
h=ufloat(6.626070040*10**(-34),0.000000081*10**(-34))
licht=299792458
lambdi=licht*h/delta

print(lambdi)
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
########################################################################################
#Ausgleichsrechnung zur SKalierung
x, y = np.genfromtxt('Messdaten/a_2_skala.txt', unpack=True)
def f(x, m, b):
    return m*x+b

params_a_2, covariance_a_2 = curve_fit(f, x, y)
errors_a_2 = np.sqrt(np.diag(covariance_a_2))
print('m_a_2 = ', params_a_2[0], '+/-', errors_a_2[0])
print('b_a_2 = ', params_a_2[1], '+/-', errors_a_2[1])
#
x_c, y_c = np.genfromtxt('Messdaten/c_skala.txt', unpack=True)
params_c, covariance_c = curve_fit(f, x_c, y_c)
errors_c = np.sqrt(np.diag(covariance_c))
print('m_c = ', params_c[0], '+/-', errors_c[0])
print('b_c = ', params_c[1], '+/-', errors_c[1])


########################################################################################
#Aufgabenteil c: Ionisierung
U_B_c, I_A_c = np.genfromtxt('Messdaten/aufg_c.txt', unpack=True)
params_ioc, covariance_ioc = curve_fit(f,params_c[0]*U_B_c, I_A_c)
errors_ioc = np.sqrt(np.diag(covariance_ioc))
print('m_ion= ', params_ioc[0], '+/-', errors_ioc[0])
print('b_ion= ', params_ioc[1], '+/-', errors_ioc[1])
M = unp.uarray(params_ioc[0], errors_ioc[0])
B = unp.uarray(params_ioc[1], errors_ioc[1])
print(-B/M)

#######################################################################################
#Differentielle Energie T = 150°C
x_, diff = np.genfromtxt('Messdaten/a_2_diff.txt', unpack=True)
U_A_ = x_ * params_a_2[0]
ascii.write([U_A_, diff],'Messdaten/franck2.tex',format='latex')

plt.plot(U_A_, diff, 'rx', label='differentiele Energieverteilung')
plt.xlabel(r'$U_{\mathrm{A}}/\si{\volt}$')
plt.ylabel(r'prop. zur Steigung')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/aufgabea2.pdf')
plt.clf()
