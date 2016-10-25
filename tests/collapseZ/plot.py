import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.22")

#set xlabel
k.xlabel = 'density/cm$^{-3}$'
k.ylabel = "Tgas/K"

k.plog = "loglog"

#plot
k.plot(columns=["rho","Tgas"], outputFileName="plot.png", blockNameLabel="Z=$Z$")

