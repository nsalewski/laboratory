import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp


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
n,f=np.genfromtxt("Messdaten/a.txt",unpack=True)
B=np.genfromtxt("Messdaten/a_B.txt",unpack=True)
G=[3,3,3,3,3]

b=f-n
g=n-35.4
err=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
error=[0.1,0.1,0.1,0.1,0.1]
b=unp.uarray(b,err)
g=unp.uarray(g,err)
bk=b[0:5]
gk=g[0:5]
B=unp.uarray(B,error)
G=unp.uarray(G,error)

V_klein=bk/gk
print(V_klein)
V_groß=B/G
ascii.write([bk,gk,B,G,V_klein,V_groß, unp.sqrt((V_klein-V_groß)**2)],"Messdaten/tab_a_B.tex",format="latex",names=["b","g","B","G","b/g","B/G","\Delta V"])

c=1/b +1/g
f_l=1/c
f_linsengleichung=ufloat(np.mean(unp.nominal_values(f_l)), np.std(unp.nominal_values(f_l)))
print("f berechnet mit Linsengleichung=", f_linsengleichung)


ascii.write([b,g, f_l],"Messdaten/tab_a.tex",format="latex",names=["b","g","f"])
for i in range(10):
    plt.plot([0,b[i].nominal_value], [g[i].nominal_value,0],label=i)


plt.xlabel(r"$g_{\mathrm{i}}$/$\si{\centi\meter}$")
plt.ylabel(r"$b_{\mathrm{i}}$/$\si{\centi\meter}$")
plt.ylim(0,20)
plt.xlim(8.5,11.0)
#plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Bilder/plot_a.pdf")
plt.clf()
a,b,c=np.genfromtxt("Messdaten/c.txt", unpack=True)
b1=unp.uarray(b-a,err)
b2=unp.uarray(b-c,err)
g1=unp.uarray(a-35.4,err)
g2=unp.uarray(c-35.4,err)
d1=(g1-b1)
d2=g2-b2
e1=g1+b1
e2=g2+b2
def f(d,e):
    return((e**2-d**2)/(4*e))
f1=f(d1,e1)
f2=f(d2,e2)
print("f1= ",f1)
print("f2= ",f2)
ascii.write([b1,g1,e1,d1,f1],"Messdaten/tab_cweiß1.tex",format="latex",names=["b1","g1","e1","d1","f1"])
ascii.write([b2,g2,e2,d2,f2],"Messdaten/tab_cweiß2.tex",format="latex",names=["b2","g2","e2","d2","f2"])

f1=ufloat(np.mean(unp.nominal_values(f1)), np.std(unp.nominal_values(f1)))
print("f1=",f1)

f2=ufloat(np.mean(unp.nominal_values(f2)), np.std(unp.nominal_values(f2)))
print("f2=",f2)
#blau
a,b,c=np.genfromtxt("Messdaten/c_blau.txt", unpack=True)
b1=unp.uarray(b-a,error)
b2=unp.uarray(b-c,error)
g1=unp.uarray(a-35.4,error)
g2=unp.uarray(c-35.4,error)
d1=(g1-b1)
d2=g2-b2
e1=g1+b1
e2=g2+b2

def f(d,e):
    return((e**2-d**2)/(4*e))
f1=f(d1,e1)
f2=f(d2,e2)
ascii.write([b1,g1,e1,d1,f1],"Messdaten/tab_cblau1.tex",format="latex",names=["b1","g1","e1","d1","f1"])
ascii.write([b2,g2,e2,d2,f2],"Messdaten/tab_cblau2.tex",format="latex",names=["b2","g2","e2","d2","f2"])

f1=ufloat(np.mean(unp.nominal_values(f1)), np.std(unp.nominal_values(f1)))
print("f1 blau=",f1)

f2=ufloat(np.mean(unp.nominal_values(f2)), np.std(unp.nominal_values(f2)))
print("f2 blau=",f2)

#rot
a,b,c=np.genfromtxt("Messdaten/c_rot.txt", unpack=True)
b1=unp.uarray(b-a,error)
b2=unp.uarray(b-c,error)
g1=unp.uarray(a-35.4,error)
g2=unp.uarray(c-35.4,error)
d1=(g1-b1)
d2=g2-b2
e1=g1+b1
e2=g2+b2

def f(d,e):
    return((e**2-d**2)/(4*e))
f1=f(d1,e1)
f2=f(d2,e2)
ascii.write([b1,g1,e1,d1,f1],"Messdaten/tab_crot1.tex",format="latex",names=["b1","g1","e1","d1","f1"])
ascii.write([b2,g2,e2,d2,f2],"Messdaten/tab_crot2.tex",format="latex",names=["b2","g2","e2","d2","f2"])

f1=ufloat(np.mean(unp.nominal_values(f1)), np.std(unp.nominal_values(f1)))
print("f1 rot =",f1)

f2=ufloat(np.mean(unp.nominal_values(f2)), np.std(unp.nominal_values(f2)))
print("f2 rot=",f2)

#############################################################################################


graw, braw, BGraw = np.genfromtxt("Messdaten/d.txt",unpack=True)
BGraw = BGraw/100
GG = 0.03
Vraw = BGraw/GG
BG = unp.uarray(BGraw, err)
V = BG / GG
g__ = graw - 35.4
b__ = braw - graw

g_ = unp.uarray(g__, err)
b_ = unp.uarray(b__, err)


ascii.write([g_, b_, V], "Messdaten/abbe.tex", format="latex")

g__ = g__/100
b__ = b__/100

def gstrich(x,f,m):
    return f*x+m

params, covariance = curve_fit(gstrich, (1+1/Vraw), g__)
errors = np.sqrt(np.diag(covariance))

plt.plot(np.linspace(0,5), params[0]*np.linspace(0,5)+params[1], 'b-',label='fit')
plt.plot((1+1/Vraw), g__, 'rx', label='Messwerte')
plt.ylabel(r"$g'$/$\si{\meter}$")
plt.xlabel(r"$(1+\frac{1}{V})$")
plt.tight_layout()
plt.savefig('Messdaten/123.pdf')

print('****************** f = ', params[0], '+/-', errors[0])
print('h = ', params[1], '+/-', errors[0])


params1, covariance1 = curve_fit(gstrich, (1+Vraw), b__)
errors1 = np.sqrt(np.diag(covariance1))

plt.clf()
plt.plot(np.linspace(0,5), (params1[0]*np.linspace(0,5)+params1[1]), 'b-', label='fit')
plt.plot((1+Vraw), b__, 'rx', label='Messwerte')
plt.ylabel(r"$b'$/$\si{\meter}$")
plt.xlabel(r"$(1+V)$")
plt.tight_layout()
plt.savefig('Messdaten/1234.pdf')
plt.clf()

print('f = ', params1[0], '+/-', errors1[0])
print('h = ', params1[1], '+/-', errors1[0])


graw1, braw1 = np.genfromtxt("Messdaten/b.txt",unpack=True)
g___ = graw1 - 35.4
b___ = braw1 - graw1
g___ = g___/100
b___ = b___/100

f = ufloat(np.mean(unp.nominal_values((1/ (1/b___ + 1/g___)))), np.std(unp.nominal_values((1/ (1/b___ + 1/g___)))))

print('f= ', f)

ascii.write([g___*100, b___*100], "Messdaten/wasser.tex", format="latex")

#########################################################################################
g___ = g___*100
b___ = b___*100

for i in range(10):
    plt.plot([0,b___[i]], [g___[i],0],label=i)

f1_=unp.uarray(params[0],errors[0])
f2_=unp.uarray(params1[0],errors1[0])
FGES = 100/f1_ + 100/f2_ - 600/(f1_*f2_)
print('FGES = ', FGES)


plt.xlabel(r"$g_{\mathrm{i}}$/$\si{\centi\meter}$")
plt.ylabel(r"$b_{\mathrm{i}}$/$\si{\centi\meter}$")
plt.ylim(0,20)
plt.xlim(5,10)
#plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Bilder/schnitti.pdf")
