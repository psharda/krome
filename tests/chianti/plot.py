import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.55")

#set xlabel
k.xlabel = "T$_{gas}$/K"
k.ylabel = "cooling/(erg cm$^{3}$ s$^{-1}$)"

k.plog = "loglog"

k.title = "KROME TEST: CHIANTI"

k.ymin = 1e-21
k.plot(columns=["Tgas","cooling"], outputFileName="plot.png",blockNameLabel="J21=$J21$")

