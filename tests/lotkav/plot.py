import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = 'time/arbitary'
k.ylabel = "population"

#plot
k.yrange([0e0,4e1])
k.plot(columns=["time","predators", "preys"],outputFileName="plot.png")

