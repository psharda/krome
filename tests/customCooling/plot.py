import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.44")

#set xlabel
k.xlabel = "T$_{gas}$/K"
#k.key = False
#k.grid = True

k.plog = "loglog"

k.title = "KROME TEST: custom cooling"

#plot
k.ylabel = "cooling/(erg cm$^{-3}$ s$^{-1}$)"

k.plot(columns=["Tgas","CUSTOM"], outputFileName="plot.png")

