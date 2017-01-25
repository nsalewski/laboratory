import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
#########
# b)
entartung, a, diff = np.genfromtxt("Messdaten/b.txt", unpack="True")
diff = diff[1:6]
diff = diff * 2  # in mm
diff = diff * 10**(-3)
print("Diff", diff)
wellenlaenge = ufloat(np.mean(diff), np.std(diff, ddof=1) / np.sqrt(len(diff)))
cadj, vnull = np.genfromtxt(
    "Messdaten/unsereMessdaten_v_null.txt", unpack="True")
v_null = ufloat(np.mean(vnull), np.std(vnull, ddof=1) / np.sqrt(len(vnull)))
c = wellenlaenge * v_null
print("Schallgeschwindigkeit=", c)
v_null_div_c = v_null / c
eins_div_lambda = 1 / wellenlaenge

######################################
# c)
delta_v = np.genfromtxt("Messdaten/adrian_clemens_e.txt", unpack="True")
# import data e, and calculate mean for each v
g1vor = ufloat(np.mean(delta_v[0:5]), np.std(
    delta_v[0:5], ddof=1) / np.sqrt(len(delta_v[0:5])))
g1rueck = ufloat(np.mean(delta_v[5:10]), np.std(
    delta_v[5:10], ddof=1) / np.sqrt(len(delta_v[5:10])))
g2vor = ufloat(np.mean(delta_v[10:15]), np.std(
    delta_v[10:15], ddof=1) / np.sqrt(len(delta_v[10:15])))
g2rueck = ufloat(np.mean(delta_v[15:20]), np.std(
    delta_v[15:20], ddof=1) / np.sqrt(len(delta_v[15:20])))
g3vor = ufloat(np.mean(delta_v[20:25]), np.std(
    delta_v[20:25], ddof=1) / np.sqrt(len(delta_v[20:25])))
g3rueck = ufloat(np.mean(delta_v[25:30]), np.std(
    delta_v[25:30], ddof=1) / np.sqrt(len(delta_v[25:30])))
g4vor = ufloat(np.mean(delta_v[30:35]), np.std(
    delta_v[30:35], ddof=1) / np.sqrt(len(delta_v[30:35])))
g4rueck = ufloat(np.mean(delta_v[35:40]), np.std(
    delta_v[35:40], ddof=1) / np.sqrt(len(delta_v[35:40])))
g5vor = ufloat(np.mean(delta_v[40:45]), np.std(
    delta_v[40:45], ddof=1) / np.sqrt(len(delta_v[40:45])))
g5rueck = ufloat(np.mean(delta_v[45:50]), np.std(
    delta_v[45:50], ddof=1) / np.sqrt(len(delta_v[45:50])))
g6vor = ufloat(np.mean(delta_v[50:55]), np.std(
    delta_v[50:55], ddof=1) / np.sqrt(len(delta_v[50:55])))
g6rueck = ufloat(np.mean(delta_v[55:60]), np.std(
    delta_v[55:60], ddof=1) / np.sqrt(len(delta_v[55:60])))
g7vor = ufloat(np.mean(delta_v[60:65]), np.std(
    delta_v[60:65], ddof=1) / np.sqrt(len(delta_v[60:65])))
g7rueck = ufloat(np.mean(delta_v[65:70]), np.std(
    delta_v[65:70], ddof=1) / np.sqrt(len(delta_v[65:70])))
g8vor = ufloat(np.mean(delta_v[70:75]), np.std(
    delta_v[70:75], ddof=1) / np.sqrt(len(delta_v[70:75])))
g8rueck = ufloat(np.mean(delta_v[75:80]), np.std(
    delta_v[75:80], ddof=1) / np.sqrt(len(delta_v[75:80])))
g9vor = ufloat(np.mean(delta_v[80:85]), np.std(
    delta_v[80:85], ddof=1) / np.sqrt(len(delta_v[80:85])))
g9rueck = ufloat(np.mean(delta_v[85:90]), np.std(
    delta_v[85:90], ddof=1) / np.sqrt(len(delta_v[85:90])))
g10vor = ufloat(np.mean(delta_v[90:95]), np.std(
    delta_v[90:95], ddof=1) / np.sqrt(len(delta_v[90:95])))
g10rueck = ufloat(np.mean(delta_v[95:100]), np.std(
    delta_v[95:100], ddof=1) / np.sqrt(len(delta_v[95:100])))
