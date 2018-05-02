import numpy as np
import matplotlib.pyplot as plt
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from uncertainties import ufloat
from astropy.io import ascii
from scipy.optimize import curve_fit

#Konstanten un Co
h = 6.626070040 *10**(-34)#Plancksche Wirkungsquantum
c = 299792458#Lichtgeschw
mu_B = 927.4009994 * 10**(-26) #Bohrsche Magneton
rot_lambda_D = 48.913 * 10**(-12) #meter
blau_lambda_D = 26.952 * 10**(-12) #meter
B_rot = 638 * 10**(-3) #milli Tesla bei 10.5 Ampere aus eichung
B_blau = 314 * 10**(-3) #milli Tesla bei 5 Ampere
#B_blau = 430 * 10**(-3) #Tesla bei 7 Ampere
B_blau_pi = 1007 *10**(-3)#Tesla bei 18 Ampere
print('##########################')
print('????????????????????????????')
print('Bei welchem Strom haben wir das Bild gemessen?')
print('5A steht doch im Dateinamen?')
print('es waren doch 1007 --milli--Tesla!?')
print('????????????????????????????')
#Rot:
print('ROT:')
lambd = 643.8 *10**(-9) #meter
r, r_links, r_rechts  = np.genfromtxt('daten/rot.txt', unpack = True)
i = 0
Delta_s = []
while i < len(r)-1:
  Delta_s = np.append(Delta_s, [r[i+1]-r[i]])
  i = i+1
print('Delta_s:')
print(Delta_s)
delta_s = []
i = 0
while i < len(r_links)-1:
    delta_s = np.append(delta_s, [r_rechts[i]-r_links[i]])
    i = i+1
print('delta_s:')
print(delta_s)
#ascii.write([Delta_s,delta_s], 'hilfsdateien/Delta_s_rot.dat', format='latex', overwrite = True)
Delta_s_mittel = ufloat(np.mean(Delta_s),np.std(Delta_s,ddof=1)/np.sqrt(len(Delta_s)))
delta_s_mittel = ufloat(np.mean(delta_s),np.std(delta_s,ddof=1)/np.sqrt(len(delta_s)))
print('Delta_s Mittelwert: ', Delta_s_mittel)
print('delta_s Mittelwert: ', delta_s_mittel)
delta_lambd = 0.5 * delta_s/Delta_s * rot_lambda_D
print('delta_lambda: ', delta_lambd)
g = (h*c* delta_lambd) / (lambd**2*mu_B*B_rot)
ascii.write([Delta_s,delta_s,delta_lambd,g], 'hilfsdateien/Delta_s_rot.dat', format='latex', overwrite = True)
g_mittel = ufloat(np.mean(g),np.std(g,ddof=1)/np.sqrt(len(g)))
print('g_mittel: ', g_mittel)
g_theo = c*h/(lambd**2)*1/4*rot_lambda_D*1/mu_B*1/B_rot
print('g_theo: ',g_theo)

#Blau:
print('######################')
print('BLAU:')
lambd = 480 * 10**(-9) #meter
#Berechnung der Abstände
b, b_links, b_rechts  = np.genfromtxt('daten/blau.txt', unpack = True)
i = 0
Delta_s = []
while i < len(b)-1:
  Delta_s = np.append(Delta_s, [b[i+1]-b[i]])
  i = i+1
print('Delta_s:')
print(Delta_s)
delta_s = []
i = 0
while i < len(b_links)-1:
    delta_s = np.append(delta_s, [b_rechts[i]-b_links[i]])
    i = i+1
print('delta_s:')
print(delta_s)
#ascii.write([Delta_s,delta_s], 'hilfsdateien/Delta_s_blau.dat', format='latex', overwrite = True)
Delta_s_mittel = ufloat(np.mean(Delta_s),np.std(Delta_s,ddof=1)/np.sqrt(len(Delta_s)))
delta_s_mittel = ufloat(np.mean(delta_s),np.std(delta_s,ddof=1)/np.sqrt(len(delta_s)))
print('Delta_s Mittelwert: ', Delta_s_mittel)
print('delta_s Mittelwert: ', delta_s_mittel)
delta_lambd = 0.5 * delta_s/Delta_s * blau_lambda_D
print('delta_lambda: ', delta_lambd)
g = (h*c* delta_lambd) / (lambd**2*mu_B*B_blau)
ascii.write([Delta_s,delta_s,np.round(delta_lambd*10**12,2),np.round(g,2)], 'hilfsdateien/Delta_s_blau.dat', format='latex', overwrite = True)
g_mittel = ufloat(np.mean(g),np.std(g,ddof=1)/np.sqrt(len(g)))
print('g_mittel: ', g_mittel)
g_theo = (c*h*blau_lambda_D)/(lambd**2*4*mu_B*B_blau)
print('g_theo: ',g_theo)

#Blau Pi
print('######################')
print('BLAU PI:')
b, b_links, b_rechts  = np.genfromtxt('daten/blau_pi.txt', unpack = True)
i = 0
Delta_s = []
while i < len(b)-1:
  Delta_s = np.append(Delta_s, [b[i+1]-b[i]])
  i = i+1
print('Delta_s:')
print(Delta_s)
delta_s = []
i = 0
while i < len(b_links)-1:
    delta_s = np.append(delta_s, [b_rechts[i]-b_links[i]])
    i = i+1
print('delta_s:')
print(delta_s)
#ascii.write([Delta_s,delta_s], 'hilfsdateien/Delta_s_blauPi.dat', format='latex', overwrite = True)
Delta_s_mittel = ufloat(np.mean(Delta_s),np.std(Delta_s,ddof=1)/np.sqrt(len(Delta_s)))
delta_s_mittel = ufloat(np.mean(delta_s),np.std(delta_s,ddof=1)/np.sqrt(len(delta_s)))
print('Delta_s Mittelwert: ', Delta_s_mittel)
print('delta_s Mittelwert: ', delta_s_mittel)
delta_lambd = 0.5 * delta_s/Delta_s * blau_lambda_D
print('delta_lambda: ', delta_lambd)
g = (h*c* delta_lambd) / (lambd**2*mu_B*B_blau_pi)
delta_lambd=np.round(delta_lambd,14)
ascii.write([Delta_s,delta_s,delta_lambd,np.round(g,2)], 'hilfsdateien/Delta_s_blauPi.dat', format='latex', overwrite = True)
g_mittel = ufloat(np.mean(g),np.std(g,ddof=1)/np.sqrt(len(g)))
print('g_mittel: ', g_mittel)
g_theo = c*h/(lambd**2)*1/4*blau_lambda_D*1/mu_B*1/B_blau_pi
print('g_theo: ',g_theo)
