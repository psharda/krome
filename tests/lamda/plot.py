import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.44")

#set xlabel
k.xlabel = "log(ntot/cm$^{-3}$)"
k.ylabel = "log(T$_{gas}$/K)"

k.plog = "loglog"
k.zplog = "log"

k.title = "KROME TEST: LAMDA"

k.plotContour(columns=["ntot","Tgas","cooling"], outputFileName="plot.png")

