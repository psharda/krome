import numpy as np
import sys
import os
import matplotlib.pyplot as plt

data_krome   = np.genfromtxt("fort.66",delimiter='  ')
data_semenov = np.genfromtxt("semenov.dat",delimiter='  ')

#Krome out
time  = np.array(data_krome[:,0])
CO    = np.array(data_krome[:,26])
HCOj  = np.array(data_krome[:,418])
H2O   = np.array(data_krome[:,58])
H2S   = np.array(data_krome[:,59])

#Semenov out
time_s  = np.array(data_semenov[:,0])
CO_s    = np.array(data_semenov[:,6])
HCOj_s  = np.array(data_semenov[:,1])
H2O_s   = np.array(data_semenov[:,2])
H2S_s   = np.array(data_semenov[:,3])

#plot
fig, ax = plt.subplots(2,2) 
#CO
ax[0,0].loglog(time, CO, label='CO KROME',marker='o', markersize='3',linestyle = 'None')
ax[0,0].loglog(time_s, CO_s, label='CO ALCHEMIC')

#HCO+
ax3 = ax[0,1].twinx()
ax3.loglog(time, HCOj, label='HCO+',marker='o',markersize='3', linestyle = 'None')
ax3.loglog(time_s, HCOj_s, label='HCO+')

ax[1,0].set_xlabel('Time [yr]')
ax[0,0].set_ylabel('fractional abundances')
ax[1,0].set_ylabel('fractional abundances')

#H2O
ax[1,0].loglog(time, H2O, label='H2O',marker='o',markersize='3', linestyle = 'None')
ax[1,0].loglog(time_s, H2O_s, label='H2O')

ax1 = ax[1,1].twinx()
ax1.loglog(time, H2S, label='H2S',marker='o',markersize='3', linestyle = 'None')
ax1.loglog(time_s, H2S_s, label='H2S')

ax[0,0].set_xlim(1e2,1.0e8)
ax[1,0].set_xlim(1e2,1.0e8)
ax[1,0].set_ylim(1e-11,5e-6)
ax[0,0].set_ylim(1e-9,1e-4)
ax[0,0].legend()

ax[0,1].set_yticklabels([])
ax3.set_xlim(1e2,1.0e8)
ax3.set_ylim(1e-13,1e-8)
ax[1,1].set_yticklabels([])
ax1.set_xlim(1e2,1.0e8)
ax1.set_ylim(1e-17,1e-10)
ax[1,1].set_xlabel('Time [yr]')
ax3.legend(loc='best')
ax1.legend(loc='best')
ax[1,0].legend(loc='best')
#plt.legend(loc='best')
#plt.savefig('plot.png', dpi=300)
plt.show()
