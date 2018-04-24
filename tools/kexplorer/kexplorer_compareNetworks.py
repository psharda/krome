import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

minAbsAbundance = 1e-10
maxAbsAbundance = 1
minRelAbundance = 1e-3
maxRelAbundance = 1e3

def compareAbundances(networkFull, networkReduced, elemInt=None, pngFolder="pngs"):

    # loop over species of interest
    if elemInt:
        for key in elemInt:
            relativeAbundanceColormap(networkFull, networkReduced, key, pngFolder)
    # loop over all species of reduced network
    else:
        for key in networkReduced.elements:
            relativeAbundanceColormap(networkFull, networkReduced, key, pngFolder)

def relativeAbundanceColormap(networkFull, networkReduced, species, pngFolder="pngs"):

    print "Making abundance colormap of %s" %(species)
    #get variable data after last time step
    x = [i[-1] for i in networkReduced.elements[species].tgasData]
    y = [i[-1] for i in networkReduced.elements[species].xvarData]
    zFull = [i[-1] for i in networkFull.elements[species].abundanceData]
    zReduced = [i[-1] for i in networkReduced.elements[species].abundanceData]

    zFull = np.asarray(zFull)
    zReduced = np.asarray(zReduced)

    zFull[zFull < minAbsAbundance] = minAbsAbundance
    zReduced[zReduced < minAbsAbundance] = minAbsAbundance

    #get relative difference
    z = np.abs(zFull - zReduced) / zFull
    # replace all extremely small numbers with a 'okay' small number
    # to avoid NaNs
    z[z < 1e-60] = 1e-60

    z = z.tolist()

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
        # logaritmic colorbar
        # Two options are presented here, and the user can choose
        # by commenting one or the other

        ### Color plot with a gradient color bar
        # rasterized=True is needed for pdfs.

        # plt.pcolormesh(x, y, z, cmap='viridis', rasterized=True,
        # norm=colors.LogNorm(vmin=zMin, vmax=zMax))

        ### Color plot with filled contours and (connected) contour lines
        # NOTE: that there is a bug in matplotlib when using 'extend':
        # "ValueError: extend kwarg does not work yet with log scale"
        # Note that there is no 'extend' arrow on the color bar...
        #
        ## OPTION 1: do a quick and dirty fix by replacing all
        # values below lower limit, with lower limit.
        #
        ## OPTION 2: wait till my pull request 10705 gets accepted
        # (https://github.com/matplotlib/matplotlib/pull/10705)
        # OR fix this in your own matplotlib version.
        # e.g. path where the file typically resides:
        # /anaconda3/envs/py27/lib/python2.7/site-packages/matplotlib/contour.py
        # Step 1: in def _init_, remove the ValueError if statement
        # Step 2: in def _process_levels,
        # change :
        # if self.extend in ('both', 'min'):
        # 	self.layers[0] = -1e150
        # to:
        # if self.extend in ('both', 'min'):
        # 	if self.logscale:
        #     self.layers[0] = 1e-150
        # 	else:
        #     self.layers[0] = -1e150
        # Reason: very large negative number gives NaN when using log,
        # so replace with very small number.
        #
        ## OPTION 1:

        # exp_min = np.floor(np.log10(zMin)-1)
        # exp_max = np.ceil(np.log10(zMax)+1)
        # lev_exp = np.arange(exp_min, exp_max)
        # levs = np.power(10, lev_exp)
        # #replace all values below lower limit, with lower limit.
        # z[z < levs[0] ] = levs[0]
        #
        # plt.contourf(x, y, z, levels=levs,
        # 			cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

        ## OPTION 2:
        exp_min = np.floor(np.log10(zMin))
        exp_max = np.ceil(np.log10(zMax)+1)
        lev_exp = np.arange(exp_min, exp_max)
        levs = np.power(10, lev_exp)
        # Try best/newest matplotlib option (OPTION 2).
        # If not present, then do OPTION 1.
        try:
            plt.contourf(x, y, z, levels=levs, extend='both',
                        cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

        except(ValueError,), err:
            # replace this with the pcolormesh option if desired.
            exp_min = np.floor(np.log10(zMin)-1)
            exp_max = np.ceil(np.log10(zMax)+1)
            lev_exp = np.arange(exp_min, exp_max)
            levs = np.power(10, lev_exp)
            #replace all values below lower limit, with lower limit.
            z[z < levs[0] ] = levs[0]

            plt.contourf(x, y, z, levels=levs,
            			cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

    else:
        #linear colorbar
        # same as in the logaritmic case:

        # plt.pcolormesh(x, y, z, cmap='viridis', rasterized=True)

        plt.contourf(x, y, z, cmap='viridis')

    #make plot labels
    plt.colorbar(label='Relative abundance', extend='both')
    plt.yscale('log')
    # plt.title('Relative abundance of %s' %(species))
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'%s (%s)' %(networkReduced.xvarName,networkReduced.xvarUnits))
    #dump png file
    print "Dumping colormap of %s" %(species)
    plt.savefig(pngFolder + '/%s' %(species))
    #plt.show()
    plt.close()
