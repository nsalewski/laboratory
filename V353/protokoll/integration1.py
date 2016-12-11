import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
x1 = np.linspace(-8, 7, 16)
x = np.linspace(-7, 7, 15)
y = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1]
z = [-0.7, 0.7, -0.7, 0.7, -0.7, 0.7, -0.7, 0.7, -
     0.7, 0.7, -0.7, 0.7, -0.7, 0.7, -0.7, 0.7]
plt.ylim(-2.5, 2.5)
plt.xlim(-6.5, 6.5)
plt.step(x, y, 'r-', label="Rechtecksfunktion")
plt.plot(x1, z, 'b-', label='Dreiecksfunktion')

plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("build/integration1.pdf")
