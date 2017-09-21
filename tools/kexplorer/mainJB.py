#this tool loads data from a file generated with
# krome_explore_flux subroutine and produces a series
# of png files to show the system evolution while evolving
import kexplorer_network

fileName = "Fluxes_demoJB.dat"
fileNameEvolution = "Grid_demoJB.dat"
oldNetwork = "react_primordial3"
newNetwork = "react_primordial3_reduced"

TgasInterest = [1e3,5e3,1e4]
xvarInterest = [1e-11,1e-13]
elemInterest = ["H","H2"]

#create network from explorer file
network = kexplorer_network.network(fileName,fileNameEvolution)

network.xvarName = "Density"
network.xvarUnits = "g/cm^3"

#find best matching criteria
#network.findBest()

#print most fluxy reactions
#network.listBest()

#dump to png folder
#network.dumpBest("pngs/")

#all below added by Jels Boulangier 30/03/2017
#create new network file for most fluxy reactions
#network.networkBest(oldNetwork,newNetwork)

#make abundace colormaps for model grid
#network.abundanceColormapAll(elemInt=elemInterest,pngFolder=pngsOut)
#if elemInt not specified, this will do all elements
#network.abundanceColormapAll(pngFolder="pngs/")


#make time evolution plot for interesting model grid points
#network.abundanceEvolution(atom="H2",tgasInt=TgasInterest,xvarInt=xvarInterest)
#network.abundanceEvolution(atom="H",tgasInt=TgasInterest,xvarInt=xvarInterest)

#ntw2tex = "/home/jels/krome/networks/react_galaxy_ism"
ntw2tex = "/media/jels/Elements/IvS/KROME/networks/Reduced_TimeDep/UMIST_manEx_TDR_1e-7Bis_noHejj.ntw"
network.network2latex(ntw2tex)
