#!/usr/bin/python
#!python


#utility to plot rate coefficients (TG, May 8 2016)
from math import log10,exp,log,sqrt
import matplotlib.pyplot as plt
import os

#reaction filename
fname = "../networks/react_COthin_rt"
#output foder
outFolder = "checkPlots"
#plot min/max temperature
TminDefault = 1e0
TmaxDefault = 1e8
#max order of magnitude to span in the plot from max value
maxSpan = 6
#number of Tgas points
imax = 100
###########################
# DO NOT EDIT BELOW
###########################


#create output folder if not present
if(not(os.path.exists(outFolder))):
	os.makedirs(outFolder)
#remove all file in folder
for filePNG in os.listdir(outFolder):
	if(filePNG=="."): continue
	if(filePNG==".."): continue
	os.unlink(outFolder+"/"+filePNG)


#expressions to be replace in Tmin/Tmax
replace = {"d":"e",".LE.":"",".GE.":"",">":"","<":""}

#shortcuts for temperature (also searched into @var)
shortcuts = {"invT":"1d0/Tgas", \
	"T32":"Tgas/3d2",\
	"T":"Tgas",\
	"invTe":"1e0/Te",\
	"sqrTgas":"sqrt(Tgas)",\
	"invTgas":"1e0/Tgas"}

#operators to replace shortcuts, e.g. *Tgas)
maths = ["+","-","*","/","(",")"]

network = dict()

#open file to read
fh = open(fname,"rb")
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	if(srow.startswith("#")): continue
	#read format
	if(srow.startswith("@format:")):
		Tmin = TminDefault
		Tmax = TmaxDefault
		idxTmin = idxTmax = idxRate = -1
		arow = [x.lower() for x in srow.split(",")]
		if("tmin" in arow): idxTmin = arow.index("tmin")
		if("tmax" in arow): idxTmax = arow.index("tmax")
		if("rate" in arow): idxRate = arow.index("rate")
		listR = [i for i in range(len(arow)) if arow[i]=="r"]
		listP = [i for i in range(len(arow)) if arow[i]=="p"]
	#store additional shortcuts (@var)
	if(srow.startswith("@var:")):
		arow = srow.split("=")
		shortcuts[arow[0].replace("@var:","").strip()] = arow[1].strip()
	if(srow.startswith("@")): continue
	arow = [x.strip() for x in srow.split(",")]
	#if Tmin is present
	if(idxTmin>-1):
		aTmin = arow[idxTmin]
		#replace symbols (e.g. <, .le.)
		for (k,v) in replace.iteritems():
			aTmin = aTmin.replace(k,v)
		#if NONE set to default min
		if(aTmin.upper()=="NONE"):
			Tmin = TminDefault
		else:
			Tmin = float(aTmin)
	#if Tmax is present
	if(idxTmax>-1):
		aTmax = arow[idxTmax]
		#replace symbols (e.g. <, .le.)
		for (k,v) in replace.iteritems():
			aTmax = aTmax.replace(k,v)
		#if NONE set to default max
		if(aTmax.upper()=="NONE"):
			Tmax = TmaxDefault
		else:
			Tmax = float(aTmax)
	#create verbatim reaction
	verbatimReaction = (" + ".join([arow[i] for i in listR if arow[i]!=""]))
	verbatimReaction += " -> "
	verbatimReaction += (" + ".join([arow[i] for i in listP if arow[i]!=""]))
	verbatimReaction = verbatimReaction.strip()
	#copy rate
	rate = arow[idxRate].lower()
	#loop for recursive shortcut replace
	for i in range(3):
		#loop on shortucts
		for k,v in shortcuts.iteritems():
			#left operator
			for mL in maths:
				#check if rate starts with shortcut+operator
				if(rate.startswith(k.lower()+mL)):
					rate = rate.replace(k.lower()+mL,"("+v.lower()+")"+mL)
				#check if rate ends with operator+shortcut
				if(rate.endswith(mL+k.lower())):
					rate = rate.replace(mL+k.lower(),mL+"("+v.lower()+")")
				#add right operator and replace if necessary
				for mR in maths:
					rate = rate.replace(mL+k.lower()+mR,mL+"("+v.lower()+")"+mR)
	#F90 -> Python math
	rate = rate.replace("dexp","exp")
	rate = rate.replace("d","e")

	#for each reaction create a dictionary with different rates and limits
	if(not(verbatimReaction in network)): network[verbatimReaction] = []
	network[verbatimReaction].append({"rate":rate,"Tmin":Tmin,"Tmax":Tmax})


irate = 0 #rate counter
outFile = outFolder+"/rateRef.dat"
fout = open(outFile,"w")
#loop on rates
for (verbatimReaction,reactions) in network.iteritems():
	plt.clf()
	nothingToPlot = True
	minKmax = minKmin = 1e99
	maxKmax = 0e0
	#loop on reactions
	for reaction in reactions:
		rate = reaction["rate"]
		Tmin = reaction["Tmin"]
		Tmax = reaction["Tmax"]
		lTmin = log10(TminDefault)
		lTmax = log10(TmaxDefault)
		xdata = []
		ydata = []
		ydata2 = []
		ydataDef = []
		#loop on Tgas points
		for ii in range(imax):
			Tgas = 1e1**(ii*(lTmax-lTmin)/(imax-1)+lTmin)
			xdata.append(Tgas)
			#try to evaluate rate for every Tgas
			try:
				kk = eval(rate.replace("tgas",str(Tgas)))
			except:
				kk = 0e0
			ydata.append(kk)
			#evaluate rate only inside limits
			if(Tgas<Tmin or Tgas>Tmax):
				ydata2.append(0e0)
			else:
				ydata2.append(kk)
				ydataDef.append(kk)
		#if rate is never evaluated rise error
		if(max(ydata)<=0e0):
			print "**********"
			print verbatimReaction
			print rate
			print "WARNING: no eval!"
			continue
		#check for negative values
		if(min(ydata)<0e0):
			print "**********"
			print rate
			print "WARNING: negative!"
		plt.loglog(xdata,ydata,"r--")
		plt.loglog(xdata,ydata2,"b")
		#store absolute min and max
		minKmax = min(minKmax,max(ydataDef))
		minKmin = min(minKmin,min(ydataDef))
		maxKmax = max(maxKmax,max(ydataDef))
		#plot junction points
		kmin = eval(rate.replace("tgas",str(Tmin)))
		kmax = eval(rate.replace("tgas",str(Tmax)))
		plt.loglog([Tmin,Tmax],[kmin,kmax],"ro")
		nothingToPlot = False

	#if there is nothing to plot skip
	if(nothingToPlot): continue
	#set default y limit
	plt.ylim(minKmin/1e1,maxKmax*1e1)
	#avoid very small values in plot
	if(minKmax/(minKmin+1e-300)>(1e1**maxSpan)):
		plt.ylim(minKmax/(1e1**maxSpan),min(maxKmax*1e1,1e0))
	#increase y axis for constant values
	if(minKmin==maxKmax):
		plt.ylim(maxKmax/1e1,maxKmax*1e1)
	plt.grid(True)
	plt.xlabel('Tgas/K')
	plt.ylabel('rate')
	plt.title(verbatimReaction+" ["+str(len(reactions))+" rate(s)]")
	#write PNG 
	fnamePNG = outFolder+"/rate"+str(int(1e6+irate))+".png"
	plt.savefig(fnamePNG)
	fout.write(fnamePNG+"\t"+verbatimReaction+"\n")
	irate += 1
print "****************"
print "****************"
print "****************"
print "output in folder "+outFolder+"/"
print "reaction list in "+outFile
print "DONE!"
fout.close()
