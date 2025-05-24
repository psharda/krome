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

Year = 3.15576e7 #seconds in a year

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

    ntot = f[1] #Get total number density since tabulated value scaled by this

    return nHnuclei

#First split the output files into separate files based on metallicity
split_by_metallicity("fort.22", "OneZoneAB")
split_by_metallicity("fort.31", "OneZoneCools")
split_by_metallicity("fort.911", "OneZoneHeats")

fig, axs = plt.subplots(nrows=3,figsize=(5.4, 14.4), tight_layout=True,sharex=True)
output_prefix = "OneZoneAB"

for i, Zval in enumerate(Zvals):
    filename = "{}_Z{}.dat".format(output_prefix, Zval)
    f = np.loadtxt(filename).T
    ab = open(filename,'r')
    nHnuclei = get_Hnuclei(filename)
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    filters = ["H","H+","H2"]
    linestyles = ["--",":","-"]
    
    cols = int(np.shape(f)[0])
    index = 0
    axs[0].plot(f[0]/(1.e6*Year), f[3], label=Zlabels[i], color=cm.lavender(i/len(Zvals)), linewidth=3.0,alpha=0.6)
    for j in range(2,cols):
        if(headers[j] in filters):
            index = np.where(np.array(filters) == headers[j])[0][0]
            if(i==0):
                axs[1].plot(f[0]/(1.e6*Year),f[j]/nHnuclei,label=str(headers[j]),color=cm.lavender(i/len(Zvals)), 
                ls=linestyles[index],linewidth=3.0,alpha=0.6)
            else:
                axs[1].plot(f[0]/(1.e6*Year),f[j]/nHnuclei,color=cm.lavender(i/len(Zvals)), 
                ls=linestyles[index],linewidth=3.0,alpha=0.6)
            index += 1

    index = 0
    filters = ["C+","C","CO"]
    for j in range(2,cols):
        if(headers[j] in filters):
            index = np.where(np.array(filters) == headers[j])[0][0]
            if(i==0):
                axs[2].plot(f[0]/(1.e6*Year),f[j]/nHnuclei/(1.6e-4*10**(Zval)),label=str(headers[j]),color=cm.lavender(i/len(Zvals)), 
                ls=linestyles[index],linewidth=3.0,alpha=0.6)
            else:
                axs[2].plot(f[0]/(1.e6*Year),f[j]/nHnuclei/(1.6e-4*10**(Zval)),color=cm.lavender(i/len(Zvals)), 
                ls=linestyles[index],linewidth=3.0,alpha=0.6)
            index += 1

axs[2].set_xlabel(r'$t \, (\rm Myr)$', fontsize=20)
axs[0].set_ylabel(r'$T_{\rm gas} \ (\rm K)$', fontsize=20)
axs[1].set_ylabel(r'$x_{i}/n_{\rm H}$', fontsize=20)
axs[2].set_ylabel(r'$x_{i}/n_{\rm C,tot}$', fontsize=20)
axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[2].set_yscale('log')

axs[1].set_ylim(1.e-6,1.2e0)
axs[2].set_ylim(1.e-6,1.2e0)

#axs.set_xscale('log')

axs[0].legend(ncol=3)
axs[1].legend(ncol=3)
axs[2].legend(ncol=3)
fig.savefig("OneZone_TimeEvolution.pdf", bbox_inches='tight')

