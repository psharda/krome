import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.24",label="dynamics")
k.loadDataFromFile("fort.34",label="chemistry")

#set xlabel
k.xlabel = "radius/cm"

k.plog = "loglog"

k.title = "KROME TEST: SHOCK 1D"

k.multiplotLayout("31",sharex=True)

k.xrange([7e17,1e19])

k.ylabel = "Tgas/K"
k.multiplotAdd(columns=["radius","Tgas"],label="dynamics")

k.ylabel = "density/(g cm$^{-3}$)"
k.multiplotAdd(columns=["radius","rho"],label="dynamics")

k.ylabel = "H$_2$/cm$^{-3}$"
k.multiplotAdd(columns=["radius","H2"],label="chemistry")

k.multiplotShow(outputFileName="plot.png")

