import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
# bisher: Berechnung der Schallgeschwindigkeit über unsere Messung mit v_0 und bestimmten lambda
# zudem: Import der ganzen daten mittels Holz-Methode, arrays für die
# difference sind jeweils fertig und können nun gegen den gang (eventuell
# auch gegen die bestimmte geschwindigkeit v des Wagens) geplottet werden


# short function for generating arrays with nominalvalues
def nomvalues_array(array):
    List = list()
    for i in range(len(array)):
        List.append(array[i].nominal_value)
    array_noms = np.asarray(List)
    return array_noms

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
# calculating diff

Gang = np.linspace(1, 10, 10)


a, vnull_adrianclemens = np.genfromtxt(
    "Messdaten/adrianundclemens/clemensadrian_v_null.txt", unpack="True")
# jedesmal muss gemessene frequenz umgerechnet werden
vnull_adrianclemens = vnull_adrianclemens * 5 / 4
vnull_ac = ufloat(np.mean(vnull_adrianclemens), np.std(
    vnull_adrianclemens, ddof=1) / np.sqrt(len(vnull_adrianclemens)))
vq = np.genfromtxt(
    "Messdaten/adrianundclemens/adrian_clemens_d.txt", unpack="True")
vq = vq * 5 / 4

g_quelle1vor = ufloat(np.mean(vq[0:5]), np.std(
    vq[0:5], ddof=1) / np.sqrt(len(vq[0:5])))
g_quelle1rueck = ufloat(np.mean(vq[5:10]), np.std(
    vq[5:10], ddof=1) / np.sqrt(len(vq[5:10])))
g_quelle2vor = ufloat(np.mean(vq[10:15]), np.std(
    vq[10:15], ddof=1) / np.sqrt(len(vq[10:15])))
g_quelle2rueck = ufloat(np.mean(vq[15:20]), np.std(
    vq[15:20], ddof=1) / np.sqrt(len(vq[15:20])))
g_quelle3vor = ufloat(np.mean(vq[20:25]), np.std(
    vq[20:25], ddof=1) / np.sqrt(len(vq[20:25])))
g_quelle3rueck = ufloat(np.mean(vq[25:30]), np.std(
    vq[25:30], ddof=1) / np.sqrt(len(vq[25:30])))
g_quelle4vor = ufloat(np.mean(vq[30:35]), np.std(
    vq[30:35], ddof=1) / np.sqrt(len(vq[30:35])))
g_quelle4rueck = ufloat(np.mean(vq[35:40]), np.std(
    vq[35:40], ddof=1) / np.sqrt(len(vq[35:40])))
g_quelle5vor = ufloat(np.mean(vq[40:45]), np.std(
    vq[40:45], ddof=1) / np.sqrt(len(vq[40:45])))
g_quelle5rueck = ufloat(np.mean(vq[45:50]), np.std(
    vq[45:50], ddof=1) / np.sqrt(len(vq[45:50])))
g_quelle6vor = ufloat(np.mean(vq[50:55]), np.std(
    vq[50:55], ddof=1) / np.sqrt(len(vq[50:55])))
g_quelle6rueck = ufloat(np.mean(vq[55:60]), np.std(
    vq[55:60], ddof=1) / np.sqrt(len(vq[55:60])))
g_quelle7vor = ufloat(np.mean(vq[60:65]), np.std(
    vq[60:65], ddof=1) / np.sqrt(len(vq[60:65])))
g_quelle7rueck = ufloat(np.mean(vq[65:70]), np.std(
    vq[65:70], ddof=1) / np.sqrt(len(vq[65:70])))
g_quelle8vor = ufloat(np.mean(vq[70:75]), np.std(
    vq[70:75], ddof=1) / np.sqrt(len(vq[70:75])))
g_quelle8rueck = ufloat(np.mean(vq[75:80]), np.std(
    vq[75:80], ddof=1) / np.sqrt(len(vq[75:80])))
g_quelle9vor = ufloat(np.mean(vq[80:85]), np.std(
    vq[80:85], ddof=1) / np.sqrt(len(vq[80:85])))
g_quelle9rueck = ufloat(np.mean(vq[85:90]), np.std(
    vq[85:90], ddof=1) / np.sqrt(len(vq[85:90])))
g_quelle10vor = ufloat(np.mean(vq[90:95]), np.std(
    vq[90:95], ddof=1) / np.sqrt(len(vq[90:95])))
g_quelle10rueck = ufloat(np.mean(vq[95:100]), np.std(
    vq[95:100], ddof=1) / np.sqrt(len(vq[95:100])))


g_quelle_rueck = [g_quelle1rueck, g_quelle2rueck, g_quelle3rueck, g_quelle4rueck, g_quelle5rueck,
                  g_quelle6rueck, g_quelle7rueck, g_quelle8rueck, g_quelle9rueck, g_quelle10rueck]
g_quelle_vor = [g_quelle1vor, g_quelle2vor, g_quelle3vor, g_quelle4vor, g_quelle5vor,
                g_quelle6vor, g_quelle7vor, g_quelle8vor, g_quelle9vor, g_quelle10vor]

# for easier array-calculation
v_null_ac_array = list()
for i in range(10):
    v_null_ac_array.append(vnull_ac)
v_null_ac_array = np.asarray(v_null_ac_array)

# arrays ready to plot, using nominalvalue for plotting
vquelle_diff_vor = unp.sqrt((g_quelle_vor - v_null_ac_array)**2)
vquelle_diff_rueck = unp.sqrt((g_quelle_rueck - v_null_ac_array)**2)
# ready to plot now

#################
vquelle_diff_vor_noms = nomvalues_array(vquelle_diff_vor)
vquelle_diff_rueck_noms = nomvalues_array(vquelle_diff_rueck)
#################

# c)part 2: difference measured by using beat-method
delta_v = np.genfromtxt(
    "Messdaten/adrianundclemens/adrian_clemens_e.txt", unpack="True")
delta_v = delta_v * 5 / 4
delta_v = delta_v / 2  # wagen bewegt sich ja mit 2v bezüglich des empfängers

# import data e, and calculate mean for each
# I know, its real crap shiiit
g1rueck = ufloat(np.mean(delta_v[0:5]), np.std(
    delta_v[0:5], ddof=1) / np.sqrt(len(delta_v[0:5])))
g1vor = ufloat(np.mean(delta_v[5:10]), np.std(
    delta_v[5:10], ddof=1) / np.sqrt(len(delta_v[5:10])))
g2rueck = ufloat(np.mean(delta_v[10:15]), np.std(
    delta_v[10:15], ddof=1) / np.sqrt(len(delta_v[10:15])))
g2vor = ufloat(np.mean(delta_v[15:20]), np.std(
    delta_v[15:20], ddof=1) / np.sqrt(len(delta_v[15:20])))
g3rueck = ufloat(np.mean(delta_v[20:25]), np.std(
    delta_v[20:25], ddof=1) / np.sqrt(len(delta_v[20:25])))
g3vor = ufloat(np.mean(delta_v[25:30]), np.std(
    delta_v[25:30], ddof=1) / np.sqrt(len(delta_v[25:30])))
g4rueck = ufloat(np.mean(delta_v[30:35]), np.std(
    delta_v[30:35], ddof=1) / np.sqrt(len(delta_v[30:35])))
g4vor = ufloat(np.mean(delta_v[35:40]), np.std(
    delta_v[35:40], ddof=1) / np.sqrt(len(delta_v[35:40])))
g5rueck = ufloat(np.mean(delta_v[40:45]), np.std(
    delta_v[40:45], ddof=1) / np.sqrt(len(delta_v[40:45])))
g5vor = ufloat(np.mean(delta_v[45:50]), np.std(
    delta_v[45:50], ddof=1) / np.sqrt(len(delta_v[45:50])))
g6rueck = ufloat(np.mean(delta_v[50:55]), np.std(
    delta_v[50:55], ddof=1) / np.sqrt(len(delta_v[50:55])))
g6vor = ufloat(np.mean(delta_v[55:60]), np.std(
    delta_v[55:60], ddof=1) / np.sqrt(len(delta_v[55:60])))
g7rueck = ufloat(np.mean(delta_v[60:65]), np.std(
    delta_v[60:65], ddof=1) / np.sqrt(len(delta_v[60:65])))
g7vor = ufloat(np.mean(delta_v[65:70]), np.std(
    delta_v[65:70], ddof=1) / np.sqrt(len(delta_v[65:70])))
g8rueck = ufloat(np.mean(delta_v[70:75]), np.std(
    delta_v[70:75], ddof=1) / np.sqrt(len(delta_v[70:75])))
g8vor = ufloat(np.mean(delta_v[75:80]), np.std(
    delta_v[75:80], ddof=1) / np.sqrt(len(delta_v[75:80])))
g9rueck = ufloat(np.mean(delta_v[80:85]), np.std(
    delta_v[80:85], ddof=1) / np.sqrt(len(delta_v[80:85])))
g9vor = ufloat(np.mean(delta_v[85:90]), np.std(
    delta_v[85:90], ddof=1) / np.sqrt(len(delta_v[85:90])))
g10rueck = ufloat(np.mean(delta_v[90:95]), np.std(
    delta_v[90:95], ddof=1) / np.sqrt(len(delta_v[90:95])))
g10vor = ufloat(np.mean(delta_v[95:100]), np.std(
    delta_v[95:100], ddof=1) / np.sqrt(len(delta_v[95:100])))


# arrays ready to plot
g_rueck = [g1rueck, g2rueck, g3rueck, g4rueck, g5rueck,
           g6rueck, g7rueck, g8rueck, g9rueck, g10rueck]

g_vor = [g1vor, g2vor, g3vor, g4vor, g5vor, g6vor, g7vor, g8vor, g9vor, g10vor]
g_rueck_noms = nomvalues_array(g_rueck)
g_rueck_vor = nomvalues_array(g_vor)
