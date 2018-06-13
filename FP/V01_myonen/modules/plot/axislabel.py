#!usr/bin/env python
#coding:utf8
import matplotlib.pyplot as plt
import numpy as np

#in rcParams wird festgelegt, dass Latex verwendet wird und siunitx in der preamble steht.
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.unicode"] = True
plt.rcParams["text.latex.preamble"].append(r"\usepackage{siunitx}")

def labels():
    axes = plt.gca() #definiert die zu verwendenden Achsen (Alle Achsen)
    x_axis = axes.get_xticks() #kriegt aktuelle Achsenbeschriftung
    label_x = [r"$\num[locale={}]{{{}}}$".format("DE", item) for item in x_axis]#formatiert Achsenbeschriftung auf neue Beschriftung mit Komma
    axes.set_xticklabels(label_x)#neue Achsenticks werden gesetzt
    y_axis = axes.get_yticks()
    label_y = [r"$\num[locale={}]{{{}}}$".format("DE", item) for item in y_axis]
    axes.set_yticklabels(label_y)


if __name__ == '__main__':
    x = np.linspace(-5, 5, 100)
    plt.plot(x, x**2, 'b-', label="test")
    labels()
    plt.savefig('axislabel.pdf')
