import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

minAbsAbundance = 1e-10
maxAbsAbundance = 1
minRelAbundance = 1e-3
maxRelAbundance = 1e3

def compareAbundances(networkFull, networkReduced, pngFolder="pngs"):

    for key in networkReduced.elements:
        relativeAbundanceColormap(networkFull, networkReduced, key, pngFolder)
    #relativeAbundanceColormap(networkFull, networkReduced, 'SIO', pngFolder)

def relativeAbundanceColormap(networkFull, networkReduced, species, pngFolder="pngs"):

    print "Making abundance colormap of %s" %(species)
    #get variable data after last time step
    x = [i[-1] for i in networkReduced.elements[species].tgasData]
    y = [i[-1] for i in networkReduced.elements[species].xvarData]
    zFull = [i[-1] for i in networkFull.elements[species].abundanceData]
    zReduced = [i[-1] for i in networkReduced.elements[species].abundanceData]

    for i, val in enumerate(zFull):
        if val < minAbsAbundance:
            zFull[i] = minAbsAbundance

    for i, val in enumerate(zReduced):
        if val < minAbsAbundance:
            zReduced[i] = minAbsAbundance

    zFull = np.asarray(zFull)
    zReduced = np.asarray(zReduced)

    #get relative difference
    z = np.abs(zFull - zReduced) / zFull
    z = z.tolist()

    for i, val in enumerate(z):
        if val == 0:
            z[i] = minRelAbundance

    #create matrix for image plot
    Ncol = len(set(x))
    Nrow = len(set(y))
    z = np.reshape(z,(Nrow, Ncol))
    x = np.reshape(x,(Nrow, Ncol))
    y = np.reshape(y,(Nrow, Ncol))

    zMin = max(z.min(), minRelAbundance)
    zMax = min(z.max(), maxRelAbundance)

    zRange = zMax/zMin

    #do no make image if abundance is too small
    if zReduced.max() < minAbsAbundance:
        print "Abundance of %s is zero everywhere" %(species)
        return

    #create image
    plt.figure()
    if(zRange>10):
        #logaritmic colorbar
        pcolormesh(x, y, z, cmap='viridis', rasterized=True,
        norm=colors.LogNorm(vmin=zMin, vmax=zMax))
    else:
        #linear colorbar
        pcolormesh(x, y, z, cmap='viridis', rasterized=True)
        
    #make plot labels
    plt.colorbar(label='Relative abundance', extend='min')
    plt.yscale('log')
    plt.title('Relative abundance of %s' %(species))
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'%s ($%s$)' %(networkReduced.xvarName,networkReduced.xvarUnits))
    #dump png file
    print "Dumping colormap of %s" %(species)
    plt.savefig(pngFolder + '/%s' %(species))
    #plt.show()
    plt.close()
