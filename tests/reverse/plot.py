import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = 'time/s'
k.ylabel = "concentration/mol"

#plot type
k.plog = "loglog"
k.ymin = 1e-6

#plot
k.plot(columns=["time","N2","O","NO","N","O2"],outputFileName="plot.png")

