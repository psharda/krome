import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.24",label="cool")
k.loadDataFromFile("shock1D.dat",label="nocool")

#set xlabel
k.xlabel = 'radius/cm'

#plot type
k.plog = "semilogy"


k.multiplotLayout("21",sharex=True)

#plot
k.ylabel = "density/(g cm$^{-3}$)"
k.yrange([1e-25,1e-23])
k.multiplotAdd(columns=["radius","rho"],label="nocool",titles="w/out cooling")
k.multiplotAddOver(columns=["radius","rho"],label="cool",titles="with cooling")

k.ylabel = "T/K"
k.yrange([5.,4e5])
k.multiplotAdd(columns=["radius","Tgas"],label="nocool",titles="w/out cooling")
k.multiplotAddOver(columns=["radius","Tgas"],label="cool",titles="with cooling")

k.multiplotShow(outputFileName="plot.png")

