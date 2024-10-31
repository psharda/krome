import numpy as np
import matplotlib.pyplot as plt
import PhysicalConstantsCGS as const
import pandas as pd

#Plot the ISM Curve at different metallicities

fig, axs = plt.subplots(nrows=2,sharex=True,figsize=(5.4,9.6),tight_layout=True)

Zvals = [-6,-5,-4,-3,-2,-1,0]
labels = [r"$10^{-6} \, Z_{\odot}$",r"$10^{-5} \, Z_{\odot}$", r"$10^{-4} \, Z_{\odot}$",r"$10^{-3} \, Z_{\odot}$",r"$10^{-2} \, Z_{\odot}$",r"$10^{-1} \, Z_{\odot}$",r"$ Z_{\odot}$"]
#colors from viridis colorbar
colors = ['#440154','#3B528B','#21918C','#5DC863','#FDE725','#FEE08B','#FC8D59']

for i,Z in enumerate(Zvals):
    data = np.loadtxt('AB_Z'+str(Z)).T
    axs[0].plot(data[0],data[1],label=labels[i],lw=3.0,alpha=0.6,color=colors[i])
    axs[1].plot(data[0],data[0]*data[1],lw=3.0,alpha=0.6,color=colors[i])

axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[0].grid()
axs[1].grid()
axs[0].set_xticks([1.e-2,1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4,1.e5,1.e6])
axs[0].set_ylabel(r'$T \, (\mathrm{K})$')
axs[1].set_ylabel(r'$P/k_{\rm B} \, (\mathrm{cm}^{-3} \mathrm{K})$')
axs[1].set_xlabel(r'$n \, (\mathrm{cm}^{-3})$')
axs[0].legend()

fig.savefig('ISMEq_Z.pdf',bbox_inches='tight')


