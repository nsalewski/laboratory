#!usr/bin/env python
#coding:utf8
#from __future__ import division
#import sys
#reload(sys)
#sys.setdefaultencoding('utf-8')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from astropy.io import ascii
from uncertainties import ufloat
import uncertainties.unumpy as unp
from modules.table import textable
import scipy.constants as const
import math as math
from modules.plot import axislabel as axis
#Daten
#rf,horizontal_1,horizontal_2, peak_1,peak_2=np.genfromtxt("data/data.txt",unpack=True)
#Fehlmessungen:
#arr1=[0.4,0.75,1.4]
#arr2=[2,3,4]
#textable.latex_tab(data=[arr1,arr2],names=[r"title column 1",r"title column 2"], filename=r"example.tex",caption=r"Beautiful caption",label=r"important_label",dec_points=[2,0])

# dec_points sets precision, i.e. dec_points[0]=2 will display 2 decimal places for all values in column 1
#Ausgleichsrechnung
#params1, covariance1 = curve_fit(theorie,rf,B1)
#errors1 = np.sqrt(np.diag(covariance1))
