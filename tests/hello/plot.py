import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = "time [arbitrary]"
k.ylabel = "abundance [arbitrary]"

k.title = "KROME TEST: HELLO"

k.plot(columns=["time","FK1","FK2","FK3"], outputFileName="plot.png")

