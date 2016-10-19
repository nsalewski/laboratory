import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


x=np.genfromtxt("1durchT.txt", unpack=True)
y=np.genfromtxt("ln(p).txt",unpack=True)
z,a=np.genfromtxt("daten_bis1bar.txt", unpack=True)

with open("daten1.txt","w") as f:
    for item in range(1,len(x)):
        f.write ("{}    &    {}    &    {}    &    {}\n".format(z[item],a[item],x[item],y[item]))
