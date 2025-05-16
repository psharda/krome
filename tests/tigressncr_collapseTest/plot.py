"""
Script to produce a plot of Tgas as a function of number density based on output from POPSICLE test
Author: Piyush Sharda (Leiden) 2024: sharda@strw.leidenuniv.nl
"""

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cm
from fileZSplit import split_by_metallicity

Zvals = [-2,-1,0]
Zlabels = [r'$10^{-2} Z_{\odot}$',r'$10^{-1} Z_{\odot}$',r'$Z_{\odot}$']

def get_Hnuclei(file):
    ab = open(file,'r')
    f = np.loadtxt(file).T
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    nHnuclei = 0.0
    cols = int(np.shape(f)[0])
    for i in range(2,cols):
        if(headers[i] == "H"):
            nHnuclei += f[i]
        elif(headers[i] == "H2"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H-"):
            nHnuclei += f[i]
        elif(headers[i] == "HD"):
            nHnuclei += f[i]
        elif(headers[i] == "OH"):
            nHnuclei += f[i]
        elif(headers[i] == "CH"):
            nHnuclei += f[i]
        elif(headers[i] == "CH2"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "CH3"):
            nHnuclei += f[i] * 3.0
        elif(headers[i] == "CH4"):
            nHnuclei += f[i] * 4.0
        elif(headers[i] == "H20_total"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H+"):
            nHnuclei += f[i]
        elif(headers[i] == "H2+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "HD+"):
            nHnuclei += f[i]
        elif(headers[i] == "OH+"):
            nHnuclei += f[i]
        elif(headers[i] == "H20+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H30+"):
            nHnuclei += f[i] * 3.0
        elif(headers[i] == "CH+"):
            nHnuclei += f[i]
        elif(headers[i] == "CH2+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H3+"):
            nHnuclei += f[i] * 3.0

    ntot = f[0] #Get total number density since tabulated value scaled by this

    return nHnuclei*ntot

#First split the output files into separate files based on metallicity
split_by_metallicity("fort.22", "CollapseAB")
split_by_metallicity("fort.31", "CollapseCools")
split_by_metallicity("fort.911", "CollapseHeats")

fig, axs = plt.subplots(ncols=1)
output_prefix = "CollapseAB"

for i, Zval in enumerate(Zvals):
    filename = "{}_Z{}.dat".format(output_prefix, Zval)
    fileAB = "CollapseAB_Z{}.dat".format(Zval)
    nH = get_Hnuclei(fileAB)
    f = np.loadtxt(filename).T
    axs.plot(nH, f[2], label=Zlabels[i], color=cm.lavender(i/len(Zvals)), linewidth=3.0,alpha=0.6)

axs.set_xlabel(r'$n_{\rm H} \, (\rm cm^{-3})$', fontsize=20)
axs.set_ylabel(r'$T_{\rm gas} \ (\rm K)$', fontsize=20)
axs.set_yscale('log')
axs.set_xscale('log')

axs.legend(ncol=3)
fig.savefig("CollapseTgas.pdf", bbox_inches='tight')

