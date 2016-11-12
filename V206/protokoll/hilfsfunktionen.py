import numpy as np

from astropy.io import ascii
t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60
T1=T1+273.15
T2=T2+273.15
pb=pb+1
pa=pa+1
def f3(t,a,alpha,b,c):
        return(((a*t**alpha)/(1+b*t**alpha)+c))
params1, covariance1 = curve_fit(f3, t, T1,bounds=([-np.inf,1,-np.inf, -np.inf],[np.inf,2,np.inf, np.inf]))
params2, covariance2 = curve_fit(f3, t, T2,bounds=([-np.inf,1,-np.inf, -np.inf],[np.inf,2,np.inf, np.inf]))
ascii.write([params1,params2,errors1,errors2], 'hilfsdateien/errparams.dat', format='latex')



errors1 = np.sqrt(np.diag(covariance1))
errors2 = np.sqrt(np.diag(covariance2))
ascii.write([t, T1, T2, pb, pa, leistung], 'hilfsdateien/values.dat', format='latex')

print(errors1)
