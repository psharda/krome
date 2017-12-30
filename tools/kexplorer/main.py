#this tool loads data from a file generated with
# krome_explore_flux subroutine and produces a series
# of png files to show the system evolution while evolving
import kexplorer_network

fileName = "demo.dat"

#create network from explorer file
network = kexplorer_network.network(fileName)
#
network.xvarName = "density"
network.xvarUnits = "g/cm^3"

#find best matching criteria
network.findBest()

#print most fluxy reactions
network.listBest()

#dump to png folder
network.dumpBest("pngs/")
