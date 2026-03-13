"""
Script to produce the evolution of a one-zone model
Author: Shyam Menon (2025); CCA, Flatiron Institute
"""

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cm
from EqTime import get_teq_onezone

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

    #ntot = f[1] #Get total number density since tabulated value scaled by this

    return nHnuclei

def get_Catoms(file):
    ab = open(file,'r')
    f = np.loadtxt(file).T
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    nCtot = 0.0
    cols = int(np.shape(f)[0])
    for i in range(2,cols):
        if(headers[i] == "C"):
            nCtot += f[i]
        elif(headers[i] == "CH"):
            nCtot += f[i]
        elif(headers[i] == "CH2"):
            nCtot += f[i]
        elif(headers[i] == "CH3"):
            nCtot += f[i]
        elif(headers[i] == "CH4"):
            nCtot += f[i]
        elif(headers[i] == "CO"):
            nCtot += f[i]
        elif(headers[i] == "CO2"):
            nCtot += f[i]
        elif(headers[i] == "C+"):
            nCtot += f[i]
        elif(headers[i] == "CO+"):
            nCtot += f[i]
        elif(headers[i] == "CH+"):
            nCtot += f[i]
        elif(headers[i] == "CH2+"):
            nCtot += f[i]

    #ntot = f[0] #Get total number density since tabulated value scaled by this

    return nCtot

fig, axs = plt.subplots(nrows=3,figsize=(5.4, 14.4), tight_layout=True,sharex=True)
output_prefix = "OneZoneAB"

colors_mpl = plt.rcParams["axes.prop_cycle"].by_key()["color"]

filename = "fort.22"
f = np.loadtxt(filename).T
ab = open(filename,'r')
nHnuclei = get_Hnuclei(filename)
header = ab.readline()
headers = header.split("#")[1].split("\n")[0].split(" ")
filters = ["H","H+","H2"]
linestyles = ["--",":","-"]
    
cols = int(np.shape(f)[0])
index = 0
axs[0].plot(f[0]/(1.e6*Year), f[3], linewidth=3.0,alpha=0.6)
for j in range(2,cols):
    if(headers[j] in filters):
        index = np.where(np.array(filters) == headers[j])[0][0]
        axs[1].plot(f[0]/(1.e6*Year),f[j]/nHnuclei,label=str(headers[j]), 
            ls=linestyles[index],color=colors_mpl[index],linewidth=3.0,alpha=0.6)
        index += 1

index = 0
filters = ["C+","C","CO"]
nCtot = get_Catoms(filename)
print(nCtot/nHnuclei)
for j in range(2,cols):
    if(headers[j] in filters):
        index = np.where(np.array(filters) == headers[j])[0][0]
        axs[2].plot(f[0]/(1.e6*Year),f[j]/nCtot,label=str(headers[j]), 
            ls=linestyles[index],linewidth=3.0,alpha=0.6,color=colors_mpl[index])
        
#Eq timescales
teq_Temp, teq_H2, teq_C, teq_Cplus, teq_CO = get_teq_onezone(filename,verbose=True)
axs[0].axvline(teq_Temp/(1.e6*Year), color=colors_mpl[0], ls='--', lw=1.5, alpha=0.6)
axs[1].axvline(teq_H2/(1.e6*Year), color=colors_mpl[2], ls='--', lw=1.5, alpha=0.6)
axs[2].axvline(teq_Cplus/(1.e6*Year), color=colors_mpl[0], ls='--', lw=1.5, alpha=0.6)
axs[2].axvline(teq_CO/(1.e6*Year), color=colors_mpl[2], ls='--', lw=1.5, alpha=0.6)

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

