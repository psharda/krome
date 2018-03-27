############## Figure settings #############
# based on matplotlib 2.0.2
# some options might not be available in older versions
import matplotlib.pyplot as plt

font = {'size'   : 40}
text = {'usetex': True}
lines = {'linewidth': 8, 'markersize': 8, 'markeredgewidth': 1}
axes = {'xmargin': 0}
#markers = {'fillstyle': 'none'}
#savefig = {'dpi': 220, 'format': 'png', 'transparent': True,'bbox': 'tight'}
savefig = {'format': 'pdf', 'transparent': True, 'bbox': 'tight'}
#figure = {'figsize': (16, 27)}
figure = {'figsize': (14, 10)}
legend = {'markerscale': 1.5, 'numpoints': 1, 'fontsize': 'xx-small', 'frameon': False, 'handlelength': 2.7}
xtick = {}
ytick = {}
plt.rc('font', **font)
#plt.rc('markers',**markers)
plt.rc('text', **text)
plt.rc('lines', **lines)
plt.rc('savefig', **savefig)
plt.rc('figure', **figure)
plt.rc('legend', **legend)
plt.rc('axes', **axes)
plt.rc('xtick', **xtick)
plt.rc('ytick', **ytick)
plt.rc('xtick.major', size=10, width=1.5)
plt.rc('xtick.minor', size=5, width=1.5)
plt.rc('ytick.major', size=10, width=1.5)
plt.rc('ytick.minor', size=5, width=1.5)
