import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = "t/yr"
k.ylabel = "fraction"

k.plog = "loglog"

k.title = "KROME TEST: cloud"

k.multiplotLayout("22")

k.yrange([1e-10,1e-3])
k.multiplotAdd(columns=["time","C"])

k.yrange([1e-12,1e-7])
k.multiplotAdd(columns=["time","OH"])

k.yrange([1e-13,1e-7])
k.multiplotAdd(columns=["time","HC3N"])

k.yrange([1e-11,1e-4])
k.multiplotAdd(columns=["time","O2"])

k.multiplotShow(outputFileName="plot.png")

