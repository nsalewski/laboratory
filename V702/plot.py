import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
nulleffekt=220
tnull=900
nulleff=220/900
zeit=ufloat(2.382,0.011)
zeit=zeit*60
print("Hwz lang=",zeit)
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

t_values = np.linspace(0,3700)
plt.plot(t_values, theorie(t_values, *params), 'g-', label="Ausgleichsgerade")
plt.errorbar(time, np.log(N), yerr=[abs(np.log(N)-np.log(N-np.sqrt(N))), abs(np.log(N+np.sqrt(N))-np.log(N))], xerr=None, fmt='rx', label="Messwerte", visible=True)
plt.ylabel(r"$\log(N(t))$/$\si{\per\second}$")
plt.xlabel(r"$t$/$\si{\second}$")
plt.ylim(0,10)
plt.xlim(0,3700)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Messdaten/indium.pdf')
plt.clf()

########Auswertung Silber##############################################################
impsges=np.genfromtxt("Messdaten/silver.txt", unpack=True)
imps=impsges[:52]
t_int=8
nullsilver=nulleff*t_int
imps=imps-nullsilver
tges=np.linspace(1,len(impsges),len(impsges))
t=np.linspace(1,len(imps),len(imps))
ascii.write([tges,np.round(impsges-nullsilver,3),np.round(np.log(impsges-nullsilver),3)],'Messdaten/anhang.tex', format='latex', names=['t','n-n0', 'ln(n)'])
print(nullsilver)
#es werden werte aus dem array entfernt, welche unter 1+nullsilver liegen, sonst probleme im logarithmus

def plot(teins,tzwei):
    arrteins=imps[0:teins]
    arreins=imps[0:teins]
    lineins=t[0:teins]
    Nzwei=imps[tzwei:]
    linzwei=t[tzwei:]
    params2lin, covariance2lin = curve_fit(theorie,8*linzwei,(Nzwei))

    arrtzwei=np.log(Nzwei)
    params2, covariance2 = curve_fit(theorie,8*linzwei,(arrtzwei))
    errors2 = np.sqrt(np.diag(covariance2))
    m2=ufloat(params2[0],errors2[0])
    b2=ufloat(params2[1],errors2[1])
    print(r"Zweiter Prozess ")
    print("m",m2,r'\n b',b2)
    m_zwei=-m2
    t_zwei=np.log(2)/m_zwei.nominal_value
    print('Halbwertszeit langer Zerfall', t_zwei)
    ascii.write([8*linzwei,np.round(Nzwei,1),np.round(np.log(Nzwei),3)],'Messdaten/tab1_silver.tex',
    format='latex', names=['t','n', 'ln(n)'])
    abzug=theorie(lineins,*params2lin)
    Neins=arrteins-abzug
    arrteins=np.log(Neins)
    params1, covariance1 = curve_fit(theorie,8*lineins,(arrteins))
    errors1 = np.sqrt(np.diag(covariance1))
    m1=ufloat(params1[0],errors1[0])
    ascii.write([8*lineins,np.round(arreins,1),np.round(Neins,2),np.round(arrteins,2)],'Messdaten/tab2_silver.tex', format='latex', names=['t','n','n-nlang','ln(n)'])
    b1=ufloat(params1[1],errors1[1])
    print(r"Erster Prozess ")
    print("m \n",m1,'b\n',b1)
    m_eins=-m1
    t_eins=np.log(2)/(m_eins.nominal_value)
    print('Halbwertszeit kurzer Zerfall', t_eins)

    linspace1=np.linspace(-10,54)

    plt.plot(linspace1*8, theorie(linspace1*8, *params1), 'g-', label=r"Ausgleichgrade Zerfall mit kleinem $T_\frac{1}{2}$ ")
    a=np.log(imps)-np.log(imps-np.sqrt(imps))
    b=(np.log(imps+np.sqrt(imps))-np.log(imps))
    for x in range(0,len(t)):
        if ((np.log(imps[x])-a[x])<=0.0):
            a[x]=np.log(imps[x])
    plt.errorbar(t*8, np.log(imps), yerr=[a,b], xerr=None, fmt='rx', label="Messwerte")
    #plt.plot(t*8, np.log(imps), 'rx', label="Messwerte")
    linspace2=np.linspace(-10,54)
    plt.plot(linspace2*8, theorie(linspace2*8,*params2),'b-', label=r"Ausgleichgrade Zerfall mit großem $T_{\frac{1}{2}}$ ")
    plt.axvline(x=teins*8, ls=':', color="g")
    plt.axvline(x=tzwei*8, ls=':', color="c", label='t*')
    plt.ylabel(r"$\ln(N(t))$")
    plt.xlabel(r"$t$/$\si{\second}$")
    plt.ylim(-0.1,5.5)
    plt.xlim(-10,54*8)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Messdaten/silber.pdf')
    plt.clf()
    an=imps-(imps-np.sqrt(imps))
    bn=((imps+np.sqrt(imps))-imps)
    plt.errorbar(t*8,imps, yerr=[an,bn], xerr=None, fmt='rx', label="Messwerte")
    t_ls=np.linspace(0, 54,100)
    plt.plot(t_ls*8,np.exp(theorie(8*t_ls,*params1))+np.exp(theorie(t_ls*8,*params2)),'g-',label="Überlagerung d. beiden Zerfälle ")

    plt.ylabel(r"$N(t)$")
    plt.xlabel(r"$t$/$\si{\second}$")
    plt.xlim(0,54*8)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Messdaten/ergebnis.pdf')
    plt.clf()
plot(13,20)



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
#
