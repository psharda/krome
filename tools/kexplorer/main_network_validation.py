# This tool is used to validate reduced networks.
# It can make absolute abundance color plots in a grid space
# for a reduced and non-reduced network. This can directly be used
# to see if both network produce the same results.
# I can also make a realtive difference abundance plot, for a more
# quantative comparison.
# It can also make abundance evolution plots for each species
# in a choosen grid space.
# These functions load data from a file that contains the abundances
# of all species at multiple time steps in a choosen grid space.
# Example format: Grid_demo.dat
import kexplorer_network
from kexplorer_compareNetworks import compareAbundances

fileName = "Fluxes_demo.dat"
# "fileNameEvolution" is a .dat files with abundances of species in temperature
# -density grid. One for the original network and one for the reduced network.
# Example file = Grid_demo.dat
fileNameEvolutionFull =
fileNameEvolutionReduced =
# When you want to plot the abundance evolution of certain species
# at choosen values of temperature and e.g. density
TgasInterest = [1e3,5e3,1e4] # in K
xvarInterest = [1e-11,1e-13] # in g/cm^3
# When you want to plot abundance color maps for a certain set of species
# in a grid of temperature and e.g. density using "fileNameEvolution"
elemInterest = ["H","H2"]

#create network objects for full and reduced from explorer file
network_reduced = kexplorer_network.network(fileName,fileNameEvolutionReduced)
network_full = kexplorer_network.network(fileName,fileNameEvolutionFull)
# Possibility to only read and store elemInterest when dealing with large networks
# network_full = kexplorer_network.network(fileName,fileNameEvolutionFull,elemInterest)

# used variable besides temperature
network_reduced.xvarName = "Density"
network_reduced.xvarUnits = "g/cm$^3$"
network_full.xvarName = "Density"
network_full.xvarUnits = "g/cm$^3$"


#### Validation of new reduced network.

# Make abundace colormaps for model grid for species of interest
# NOTE: the user can choose between two kind of colormaps
# see 'abundanceColormap()' in kexplorer_network for details
network_full.abundanceColormapAll(elemInt=elemInterest, timeEvolution=True, pngFolder="tpngsFull/")
network_reduced.abundanceColormapAll(elemInt=elemInterest,pngFolder="tpngsRed/")
# If elemInt not specified, this will do all species
# network.abundanceColormapAll(pngFolder="pngs/")
# timeEvolution is an option to make colormaps for every time step to visualise
# temporal evolution in temparature-rho space

# Make a video of all time evolution color plots for species of interest
network_full.makeEvolutionVideo(elemInt=elemInterest, pngFolder="tpngsFull/")



# Make relative abundace colormaps between full and reduced network
compareAbundances(network_full, network_reduced,
                pngFolder="tpngsCompare/")
# Do this only for species in elemInt
# kexplorer_compareNetworks.compareAbundances(network_full, network_reduced,
#                                             elemInt=elemInterest,
#                                             pngFolder="tpngsCompare/")


#make time evolution plot for interesting model grid points
network_full.abundanceEvolution(species=["H2"], tgasInt=TgasInterest,
                                xvarInt=xvarInterest, pngFolder="pngsFull/")
network_full.abundanceEvolution(species=["H"], tgasInt=TgasInterest,
                                xvarInt=xvarInterest, pngFolder="pngsFull/")
network_red.abundanceEvolution(species=["H2"], tgasInt=TgasInterest,
                                xvarInt=xvarInterest, pngFolder="pngsRed/")
network_red.abundanceEvolution(species=["H"], tgasInt=TgasInterest,
                                xvarInt=xvarInterest, pngFolder="pngsRed/")
#make time evolution plot with multiple species on one plot
network_full.abundanceEvolution(species=elemInterest, tgasInt=TgasInterest,
                                xvarInt=xvarInterest, pngFolder="pngsFull/")
