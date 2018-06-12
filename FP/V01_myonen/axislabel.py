#!usr/bin/env python
#coding:utf8
import matplotlib.pyplot as plt
import numpy as np
#print(plt.rcParams["text.latex.preamble"])
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.unicode"] = True
plt.rcParams["text.latex.preamble"].append(r"\usepackage{siunitx}")
def labels():
    axes = plt.gca() #definiert die zu verwendenden achsen
    x_axis = axes.get_xticks() #kriegt aktuelle achsenbeschriftung
    label_x = [r"$\num[locale={}]{{{}}}$".format("DE", item) for item in x_axis]#formatiert achsenbeschriftung auf neue beschriftung mit komma
    axes.set_xticklabels(label_x)
    y_axis = axes.get_yticks()
    label_y = [r"$\num[locale={}]{{{}}}$".format("DE", item) for item in y_axis]
    axes.set_yticklabels(label_y)

x = np.linspace(-2, 2, 1000)
plt.plot(x, x**2, 'b-', label="test")
labels()
plt.savefig('test.pdf')
