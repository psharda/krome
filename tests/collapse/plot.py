import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.22")

#set xlabel
k.xlabel = r'density/cm$^{-3}$'

#plot type
k.plog = "loglog"

#plot
k.ylabel = "fraction"
k.ymax = 3e0
k.multiplotLayout("21",sharex=True)
k.multiplotAdd(columns=["ntot","H","H2"])

k.ymin = 1e2
k.ylabel = "T/K"
k.multiplotAdd(columns=["ntot","Tgas"])
k.multiplotShow(outputFileName="plot.png")

