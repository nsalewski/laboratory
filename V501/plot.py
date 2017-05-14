import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
def y(x, m, b):
    return m * x + b

##########################################################################################
# E-Feld
x=np.linspace(-12,38)
n_, v_säge = np.genfromtxt("Messdaten/frequenzsaege.txt",unpack=True)
ascii.write([n_, v_säge], 'Messdaten/tab_saegi.tex', format="latex", names=['Frequenzverhältnis','frequenz'])
vwechsel=v_säge/n_
vwechsel=ufloat(np.mean(vwechsel),np.std(vwechsel, ddof=1) / np.sqrt(len(vwechsel)))
print(vwechsel)
D, Ud400, Ud300, Ud200 = np.genfromtxt("Messdaten/efeld.txt",unpack=True)
ascii.write([D*2.54, Ud400, Ud300, Ud200], 'Messdaten/tab_efeld.tex', format="latex")
D=D*2.54
params400, covariance400 = curve_fit(y,Ud400,D)
errors400 = np.sqrt(np.diag(covariance400))
params300, covariance300 = curve_fit(y, Ud300,D)
errors300 = np.sqrt(np.diag(covariance300))
params200, covariance200 = curve_fit(y, Ud200,D)
errors200 = np.sqrt(np.diag(covariance200))
print('m400 = ', params400[0], '+/-', errors400[0])
print('m300 = ', params300[0], '+/-', errors300[0])
print('m200 = ', params200[0], '+/-', errors200[0])


m=[params200[0],params300[0],params400[0]]
Ud=[10**3/200,10**3/300,10**3/400]
paramsud, covarianceud = curve_fit(y,Ud,m)
errorsud = np.sqrt(np.diag(covarianceud))
print('m_ud = ', paramsud[0], '+/-', errorsud[0])
Uud=np.linspace(1/160,1/460)
Uud=Uud*10**3
plt.plot(Uud,paramsud[0]*Uud+paramsud[1], 'b-',label=r'Regressionsgrade')
plt.plot(Ud,m, 'rx', label=r'Messwerte')
plt.ylabel(r"$\frac{D}{U_\mathrm{d}}$/$\si{\centi\meter\per\volt}$")
plt.xlabel(r"$\frac{1}{U_\mathrm{B}}\cdot 10^3$/$\si{\per\volt}$")
plt.xlim(2.2,6.0)
#plt.ylim(-2,14)
plt.legend()
plt.tight_layout()
plt.savefig('Messdaten/plotm.pdf')
plt.clf()

plt.plot(x, params200[0]*x+params200[1], 'g-',label=r'Regression $U_\mathrm{B}=\SI{200}{Volt}$')
plt.plot(Ud200,D, 'gx', label=r'Messwerte $U_\mathrm{B}=\SI{200}{Volt}$')
plt.plot(x, params300[0]*x+params300[1], 'b-',label=r'Regression $U_\mathrm{B}=\SI{300}{Volt}$ ')
plt.plot(Ud300,D, 'bx', label=r'Messwerte $U_\mathrm{B}=\SI{300}{Volt}$')
plt.plot(x, params400[0]*x+params400[1], 'r-',label=r'Regression $U_\mathrm{B}=\SI{400}{Volt}$ ')
plt.plot(Ud400,D, 'rx', label=r'Messwerte $U_\mathrm{B}=\SI{400}{Volt}$')
plt.ylabel(r"$D$/$\si{\centi\meter}$")
plt.xlabel(r"$U_\mathrm{d}$/$\si{\volt}$")
plt.xlim(-12,38)
plt.ylim(-2,14)
plt.legend()
plt.tight_layout()
plt.savefig('Messdaten/plotefeld.pdf')
plt.clf()
#########################################################################################
# B-Feld

I250, D_, I450  = np.genfromtxt("Messdaten/messdaten502a.txt",unpack=True)
ascii.write([D_*2.54, I250, I450], 'Messdaten/tab_bfeld.tex', format="latex")


params, covariance = curve_fit(y, 4*np.pi*10**(-7)*8/np.sqrt(125)*20*I250/0.282, D_/(D_**2+0.143**2))
errors = np.sqrt(np.diag(covariance))

print('m = ', params[0], '+/-', errors[0])
print('b = ', params[1], '+/-', errors[1])

m = unp.uarray(params[0], errors[0])
e_theo = unp.uarray(1.6021766208*10**(-19), 0.0000000098*10**(-19))
m_theo = unp.uarray(9.10938356*10**(-31), 0.00000011*10**(-31))
e_m=m**2*8*250
e_m_theo = e_theo/m_theo
print('experiment = ', e_m)
print('theorie = ', e_m_theo)

plt.plot(np.linspace(0,0.0002), params[0]*np.linspace(0,0.0002)+params[1], 'b-',label='Ausgleichsgerade')
plt.plot(4*np.pi*10**(-7)*8/np.sqrt(125)*20*I250/0.282,D_/(D_**2+0.143**2) , 'rx', label='Messwerte')
plt.ylabel(r"$\frac{D}{D^2+L^2}$/$\si{\per\meter}$")
plt.xlabel(r"$B$/$\si{\tesla}$")
plt.tight_layout()
plt.savefig('Messdaten/plotbfeld.pdf')
plt.clf()

D_ = D_[0:-1]
I450 = I450[0:-1]
params, covariance = curve_fit(y, 4*np.pi*10**(-7)*8/np.sqrt(125)*20*I450/0.282, D_/(D_**2+0.143**2))
errors = np.sqrt(np.diag(covariance))

print('m = ', params[0], '+/-', errors[0])
print('b = ', params[1], '+/-', errors[1])

plt.plot(np.linspace(0,0.0002), params[0]*np.linspace(0,0.0002)+params[1], 'b-',label='Ausgleichsgerade')
plt.plot(4*np.pi*10**(-7)*8/np.sqrt(125)*20*I450/0.282,D_/(D_**2+0.143**2) , 'rx', label='Messwerte')
plt.ylabel(r"$\frac{D}{D^2+L^2}$/$\si{\per\meter}$")
plt.xlabel(r"$B$/$\si{\tesla}$")
plt.tight_layout()
plt.savefig('Messdaten/plotbfeld2.pdf')
















#plt.plot(theta, w/1000, 'rx', label="Messwerte")
#plt.plot(thetaplot, theorie(thetaplot)/1000, 'b-', label="Theoriekurve")
#
#plt.ylabel(r"$\omega/\si{\kilo\hertz}$")
#plt.xlabel(r"$\theta/\si{\radian}$")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig('Bilder/b1.pdf')
#
