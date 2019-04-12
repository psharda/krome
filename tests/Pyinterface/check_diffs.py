import numpy as np
import matplotlib.pyplot as plt

files = [["idoff.f.dat", "idoff.py.dat"],
         ["idon.f.dat", "idon.py.dat"]]

colors =  ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#CFECF9', '#7F7F7F', '#BCBD22', '#17BECF']
markers = [None, "x"]
linewidths = [1, 0]

for fnames in files:
    plt.clf()
    for ii, fname in enumerate(fnames):
        data = np.loadtxt(fname).T
        for jj, dd in enumerate(data[2:]):
            label = ""
            if jj == 0:
                label = fname
            plt.loglog(data[0],
                       dd,
                       color=colors[jj],
                       marker=markers[ii],
                       markersize=5.,
                       linewidth=linewidths[ii],
                       label=label)

    plt.ylim(1e-8, 3e0)
    plt.ylabel("abundances / cm-3")
    plt.xlabel("time / yr")
    plt.legend(loc="best")
    plt.show()

