import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
x = np.linspace(-10, 10, 300)
plt.plot(x, np.sin(x), 'r-', label="Sinus")
plt.plot(x, -0.7 * np.cos(x), 'b-', label='Cosinus')
plt.xlim(-10, 10)
plt.xlabel(r"$x$")
plt.ylabel(r"$f(x)$")
plt.legend(loc='upper right')
plt.grid()
plt.tight_layout()
plt.savefig("build/integration2.pdf")
