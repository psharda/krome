import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = 'time/yr'
k.ylabel = "abundance/cm$^{-3}$"

k.plog = "loglog"

#plot
k.yrange([1e-2,2e0])
k.oplot(columns=["time","C"], blockNameLabel="C, J21=$j21$")
k.plot(columns=["time","C+"], outputFileName="plot.png", blockNameLabel="C+, J21=$j21$")

