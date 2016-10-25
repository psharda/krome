import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.24")

#set xlabel
k.xlabel = "radius/cm"
#k.key = False
#k.grid = True

k.multiplotLayout("12")
k.plog = "loglog"

k.title = "KROME TEST: compact"

#plot
k.ylabel = "T$_{gas}$"
k.multiplotAdd(columns=["radius","Tgas"])

k.ylabel = "density/(g/cm$^3$)"
k.multiplotAdd(columns=["radius","density"])
k.multiplotShow(outputFileName="plot.png")

