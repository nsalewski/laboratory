import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
nulleffekt=220
tnull=900
nulleff=220/900
def theorie(x,m,b):
    return m*x+b


#####################################################################################
#Indium
number, N = np.genfromtxt('Messdaten/indium.txt', unpack=True)
time_interval = 240
time = number*240
ascii.write([time, N], 'Messdaten/tab_indium.tex', format="latex")
N = N - 176/3

params, covariance = curve_fit(theorie, time, np.log(N))
errors = np.sqrt(np.diag(covariance))
pace = ufloat(params[0],errors[0])

print('m = ', params[0], ' +/- ', errors[0], ' b = ', params[1], ' +/- ', errors[1])
print('Halbertszeit = ', np.log(2) / pace)






########Auswertung Silber##############################################################
imps=np.genfromtxt("Messdaten/silver.txt", unpack=True)
imps=imps[:52]
t_int=8
nullsilver=nulleff*t_int
imps=imps-nullsilver
t=np.linspace(0,len(imps),len(imps))
print(nullsilver)
#es werden werte aus dem array entfernt, welche unter 1+nullsilver liegen, sonst probleme im logarithmus

#def plot(teins,tzwei):
#    arrteins=imps[0:teins]
#    lineins=t[0:teins]
#    arrtzwei=imps[tzwei:len(imps)]
#    linzwei=t[tzwei:len(imps)]
#    params2lin, covariance2lin = curve_fit(theorie,8*linzwei,(arrtzwei))
#
#    arrtzwei=np.log(arrtzwei)
#    print(arrtzwei)
#    params2, covariance2 = curve_fit(theorie,8*linzwei,(arrtzwei))
#    errors2 = np.sqrt(np.diag(covariance2))
#    m2=ufloat(params2[0],errors2[0])
#    b2=ufloat(params2[1],errors2[1])
#    print(r"Zweiter Prozess ")
#    print("m",m2,r'\n b',b2)
#    m_zwei=-m2
#    t_zwei=np.log(2)/m_zwei.nominal_value
#    print('Halbwertszeit langer Zerfall', t_zwei)
#
#    abzug=theorie(lineins,*params2lin)
#    arrteins=arrteins-abzug
#    arrteins=np.log(arrteins)
#    params1, covariance1 = curve_fit(theorie,8*lineins,(arrteins))
#    errors1 = np.sqrt(np.diag(covariance1))
#    m1=ufloat(params1[0],errors1[0])
#
#    b1=ufloat(params1[1],errors1[1])
#    print(r"Erster Prozess ")
#    print("m",m1,r'\n b',b1)
#    m_eins=-m1
#    t_eins=np.log(2)/(m_eins.nominal_value)
#    print('Halbwertszeit langer Zerfall', t_eins)
#
#    linspace1=np.linspace(-10,54)
#
#    plt.plot(linspace1*8, theorie(linspace1*8, *params1), 'g-', label="Ausgleichsgrade 1")
#    a=np.log(imps)-np.log(imps-np.sqrt(imps))
#    b=(np.log(imps+np.sqrt(imps))-np.log(imps))
#    for x in range(0,len(t)):
#        if ((np.log(imps[x])-a[x])<=0.0):
#            a[x]=np.log(imps[x])
#    plt.errorbar(t*8, np.log(imps), yerr=[a,b], xerr=None, fmt='x', label="Messwerte")
#    #plt.plot(t*8, np.log(imps), 'rx', label="Messwerte")
#    linspace2=np.linspace(-10,54)
#    plt.plot(linspace2*8, theorie(linspace2*8,*params2),'b-', label="Ausgleichsgrade 2")
#    plt.axvline(x=teins*8, ls=':', color="g")
#    plt.axvline(x=tzwei*8, ls=':', color="c", label=r"t^{*}")
#    plt.ylabel(r"$\log(N(t))$/$\si{\per\second}$")
#    plt.xlabel(r"$t$/$\si{\second}$")
#    plt.ylim(-0.1,5.5)
#    plt.xlim(-10,54*8)
#    plt.legend(loc='best')
#    plt.tight_layout()
#    plt.savefig('Messdaten/silber.pdf')
#    plt.clf()
#    plt.errorbar(t*8, np.log(imps), yerr=[a,b], xerr=None, fmt='x', label="Messwerte")
#    plt.plot(t*8,np.exp(-m_eins.nominal_value*8*t)+np.exp(-m_zwei.nominal_value*8*t),'g-',label="Thei")
#    plt.ylabel(r"$\log(N(t))$/$\si{\per\second}$")
#    plt.xlabel(r"$t$/$\si{\second}$")
#    plt.legend(loc='best')
#    plt.tight_layout()
#    plt.savefig('Messdaten/ergebnis.pdf')
#    plt.clf()
#plot(13,25)
#


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
