import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.22")

#set xlabel
k.xlabel = "density/cm$^{-3}$"

k.plog = "loglog"

k.title = "KROME TEST: collapse UV"

k.multiplotLayout("21")


k.ylabel = "T/K"
k.keyPosition("br")
k.multiplotAdd(columns=["ntot","Tgas"],blockNameLabel="J21=$J21$")

k.ylabel = "xH$_2$"
k.ymin = 1e-8
k.key = False
k.multiplotAdd(columns=["ntot","H2"])

k.multiplotShow(outputFileName="plot.png")

