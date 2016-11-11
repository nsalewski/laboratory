import matplotlib.pyplot as plt
import numpy as np


t,T1,T2,pb,pa,leistung =np.genfromtxt("daten.txt", unpack=True)
t=t*60
plt.plot(t, T1, 'o', label='Messdaten T1 (warmes Reservoir)')
plt.plot(t, T2, 'o', label= 'Messdaten T2 (kaltes Reservoir')
plt.savefig('build/plot.pdf')
