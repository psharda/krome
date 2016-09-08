import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.22")

#set xlabel
k.xlabel = "density/cm$^{-3}$"
#k.key = False
#k.grid = True

k.plog = "loglog"

k.title = "KROME TEST: collapse UV+X"

#plot
k.ylabel = "T$_{gas}$/K"

k.plot(columns=["ntot","Tgas"],outputFileName="plot.png",blockNameLabel="J21=$0$")

