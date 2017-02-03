import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


g = ufloat(9.811899, 0.000041)
l = ufloat(0.720, 0.003)

T_einzeln = np.genfromtxt('Messdaten/b.txt', unpack='True')
T_einzeln = T_einzeln/5

T_links = T_einzeln[0:10]
T_rechts = T_einzeln[10:20]

T_gleichphasig = T_einzeln[20:30]
T_gegenphasig = T_einzeln[30:40]

T_schwingung, T_schwebung = np.genfromtxt('Messdaten/schwebungb.txt', unpack='True')
T_schwebung = T_schwebung/5

########################################################################################

ascii.write([T_links, T_rechts], 'Messdaten/einzeln_b.tex', format='latex')
ascii.write([T_gleichphasig, T_gegenphasig], 'Messdaten/gekoppelt_b.tex', format='latex')
ascii.write([T_schwingung, T_schwebung], 'Messdaten/schwebung_b.tex', format='latex')

########################################################################################

T_links = ufloat(np.mean(T_links), np.std(T_links, ddof=1)/np.sqrt(len(T_links)))
T_rechts = ufloat(np.mean(T_rechts), np.std(T_rechts, ddof=1)/np.sqrt(len(T_rechts)))
T_gleichphasig = ufloat(np.mean(T_gleichphasig), np.std(T_gleichphasig, ddof=1)/np.sqrt(len(T_gleichphasig)))
T_gegenphasig = ufloat(np.mean(T_gegenphasig), np.std(T_gegenphasig, ddof=1)/np.sqrt(len(T_gegenphasig)))
T_schwebung = ufloat(np.mean(T_schwebung), np.std(T_schwebung, ddof=1)/np.sqrt(len(T_schwebung)))
T_schwingung = ufloat(np.mean(T_schwingung), np.std(T_schwingung, ddof=1)/np.sqrt(len(T_schwingung)))
w_plus_theo = unp.sqrt(g/l)
w_plus_exp = 2*np.pi / T_gleichphasig
w_minus_exp = 2*np.pi / T_gegenphasig
T_plus_theo = 2*np.pi*unp.sqrt(l/g)
kappa = (T_gleichphasig**2 - T_gegenphasig**2)/(T_gleichphasig**2 + T_gegenphasig**2)
w_minus_theo = unp.sqrt(w_plus_theo**2*(1+kappa)/(1-kappa))
T_minus_theo = 2*np.pi/w_minus_theo

########################################################################################

print('T_links = ', T_links)
print('T_rechts = ', T_rechts)
print('T_gleichphasig = ', T_gleichphasig)
print('T_gegenphasig = ', T_gegenphasig)
print('w_plus_exp = ', w_plus_exp)
print('w_minus_exp = ', w_minus_exp)
print('w_plus_theo = ', w_plus_theo)
print('T_plus_theo', T_plus_theo)
print('kappa = ', kappa)
print('w_minus_theo = ',  w_minus_theo)
print('T_minus_theo = ', T_minus_theo)
print('T_schwingung = ', T_schwingung)
print('T_schwebung = ', T_schwebung)
print('T_schwebung_theo = ', T_plus_theo*T_minus_theo/(T_plus_theo-T_minus_theo))
print('w_theo_1 = ', w_plus_theo - w_minus_theo)
print('w_exp_1 = ', w_plus_exp - w_minus_exp)
print('w_theo_2 = ', 2*np.pi/(T_plus_theo*T_minus_theo/(T_plus_theo-T_minus_theo)))
print('w_exp_1 = ', 2*np.pi/T_schwebung)
