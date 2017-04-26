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
ascii.write([bk,gk,B,G,V_klein,V_groß],"Messdaten/tab_a_B.tex",format="latex",names=["b","g","B","G","b/g","B/G"])
V_klein=ufloat(np.mean(unp.nominal_values(V_klein)), np.std(unp.nominal_values(V_klein)))
V_groß=ufloat(np.mean(unp.nominal_values(V_groß)), np.std(unp.nominal_values(V_groß)))
print("V_klein=", V_klein)
print("V_groß=", V_groß)

c=1/b +1/g
f_l=1/c
f_linsengleichung=ufloat(np.mean(unp.nominal_values(f_l)), np.std(unp.nominal_values(f_l)))
print("f berechnet mit Linsengleichung=", f_linsengleichung)


ascii.write([b,g, f_l],"Messdaten/tab_a.tex",format="latex",names=["b","g","f"])
for i in range(10):
    plt.plot([0,b[i].nominal_value], [g[i].nominal_value,0],label=i)


plt.xlabel(r"$g_{\mathrm{i}}$")
plt.ylabel(r"$b_{\mathrm{i}}$")
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
