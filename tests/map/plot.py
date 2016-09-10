import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66")

#set xlabel
k.xlabel = 'log[density/(g cm$^{-3}$)]'
k.ylabel = "log(J21)"

k.plog = "loglog"
k.zplog = "log"

#plot
k.plotContour(columns=["rho","J21","H2"], outputFileName="plot.png")

