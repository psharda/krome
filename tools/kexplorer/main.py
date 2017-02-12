#this tool loads data from a file generated with
# krome_explore_flux subroutine and produces a series
# of png files to show the system evolution while evolving
import kexplorer_network

fileName = "demo.dat"

#create network from explorer file
network = kexplorer_network.network(fileName)

#find best matching criteria
network.findBest(atom="H")

#print most fluxy reactions
network.listBest()

#dump to png folder
network.dumpBest()

