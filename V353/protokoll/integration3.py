import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
z=np.linspace(-1,1)

c=[-10,-8,-6,-4,-2,0,2,4,6]
y = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1]
for item in c:
        x=np.linspace(item+1,item+2)
        plt.plot(x+1, z**2-1, "r-")
        x_2=np.linspace(item,item+1)
        plt.plot(x_2+1,-z**2+1,"r-")
x1 = np.linspace(-8,7,16)
z=[-1.6,1.6,-1.6,1.6,-1.6,1.6,-1.6,1.6,-1.6,1.6,-1.6,1.6,-1.6,1.6,-1.6,1.6]
plt.ylim(-2.5,2.5)
plt.xlim(-3,3)
plt.plot(x1,z, 'b-', label='Dreiecksfunktion')
plt.plot(13,233**2+1,"r-", label="Zusammengesetzte Funktion aus Parabeln")
plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("build/integration3.pdf")
