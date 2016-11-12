import numpy as np

from astropy.io import ascii
t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60
T1=T1+273.15
T2=T2+273.15
pb=pb+1
pa=pa+1
ascii.write([t, T1, T2, pb, pa, leistung], 'values.dat', format='latex')
