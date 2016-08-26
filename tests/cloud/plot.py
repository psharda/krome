import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = "time/yr"

#plot type
k.plog = "loglog"

#plot
k.ylabel = "fraction"
k.multiplotLayout("22")
k.multiplotAdd(columns=["time","C"])
k.multiplotAdd(columns=["time","OH"])
k.multiplotAdd(columns=["time","HC3N"])
k.multiplotAdd(columns=["time","O2"])
k.multiplotShow()

