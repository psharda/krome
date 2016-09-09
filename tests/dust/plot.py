import kp

k = kp.kp()

#load data
k.loadDataFromFile("fort.66",label="gas")
k.loadDataFromFile("fort.77",label="dust")

#set xlabel
k.xlabel = "T$_{gas}$/K"
#k.key = False
#k.grid = True

k.plog = "loglog"

k.title = "KROME TEST: dust"

#plot
k.ylabel = "$\\rho$/(g cm$^{-3}$)"
k.ymin = 1e-25

k.oplot(columns=["Tgas","mC","mSi"],label="gas")
k.oplot(columns=["Tgas","rhodust"],label="dust", plotStyle=["b-"],condition="($dustType$==0)",keyLabelAll="dust C")
k.plot(columns=["Tgas","rhodust"],label="dust", plotStyle=["r-"],condition="($dustType$==1)",keyLabelAll="dust Si",outputFileName="plot.png")

