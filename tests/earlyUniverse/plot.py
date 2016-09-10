import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.22")

#set xlabel
k.xlabel = 'log(1+z)'
k.ylabel = "n/ntot"

k.plog = "loglog"

k.multiplotLayout("22",sharex=True)
#plot
k.yrange([1e-25,1e1])
k.multiplotAdd(columns=["z","H", "H+", "H-", "H2", "H2+", "H3+"])
k.multiplotAdd(columns=["z","HE", "HE+", "HE++", "HEH+"])

k.multiplotAdd(columns=["z","D", "D+", "HD", "HD+", "H2D+"])

k.ylabel = "Tgas/K"
k.autoscale()
k.multiplotAdd(columns=["z","Tgas", "Trad"])

k.multiplotShow(outputFileName="plot.png")

