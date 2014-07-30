#!/usr/bin/python
#!python

import sys
from os.path import exists as file_exists
from math import log,log10,exp,sqrt
from scipy.interpolate import interp1d as interp_spline
from os import listdir
from os.path import isfile, join,isdir
import numpy as np
import matplotlib.pyplot as plt

foutname = "coolChianti.dat"
path = "chianti/"
plotpre = "none" #"c_2"
levLimit = 20
flist = []
dirs = [f for f in listdir(path) if isdir(join(path,f))]
usemetal = ["c","o","ne","fe"]
for dr in dirs:
	mypath = path+dr+"/"
	dirs2 = [f for f in listdir(mypath) if isdir(join(mypath,f))]
	flist += [mypath+x for x in dirs2]

#flist = ["chianti/c/c_2"]

##################################
#return boolean true if argument is number
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

##################################
#double format snum
def format_double(snum):
	snum = str(snum).lower()
	#format string number to F90 double
	if("d" in snum): return snum
	if("e" in snum): return snum.replace("e","d")	
	return snum+"d0"

######################################
def read_kex(fname,collname):
	#****SPLUPS****
	#method: Burgess+Tully 1992 A&A
	if(not(file_exists(fname))): 
		print("WARNING: no data for coolision with "+collname+", "+fname+" does not exists!")
		return
	if(collname=="e"):
		iblocks = 5 #number of integer blocks of size 3 characters
		boff = 0 #block offset, since H+ collision misses the 2 initial integers
	elif(collname=="H+"):
		iblocks = 3 #number of integer blocks of size 3 characters
		boff = 2 #block offset, since H+ collision misses the 2 initial integers
	else:
		sys.exit("ERROR: ",collname+" collider not recognized in read_kex!")
	
	print "reading "+fname
	fh = open(fname,"rb") #open file
	fout.write("\n#collider, level_up, level_down, rate\n")
	itrans = 0
	allrate = dict()
	for row in fh:
		srow = row.strip()
		if(srow==""): continue #skip blank
		if(srow[0]=="%"): break #if comment breaks
		if(srow=="-1"): break #if end data breaks
		arow = [row[j*3:(j+1)*3] for j in range(iblocks)]
		itrans += 1
		for j in range(len(row[iblocks*3:])/10):
			block = [row[iblocks*3+j*10:iblocks*3+(j+1)*10].strip()]
			if(block!=[""]): arow += block
		#arow += [x for x in row[15:].strip().split(" ") if is_number(x)] #keep only columns with numbers
		levLow = int(arow[2-boff])
		levUp = int(arow[3-boff])
		if(levLow>levLimit or levUp>levLimit): continue
		if(plotpre in pre): l = plt.axvline(x=itrans+.5, ymin=float(levLow-1.)/(maxlev+1), ymax=float(levUp-1.)/(maxlev+1),color="r")
		ty = int(arow[4-boff]) #type
		dE = float(arow[6-boff]) #Ry
		cc = float(arow[7-boff]) #scaling factor
		spli = [float(x) for x in arow[8-boff:]] #read spline
		xx = [float(x)/(len(spli)-1) for x in range(len(spli))] #normalized x value for spline
		imax = int(1e1) #number of Tgas bins
		ff = interp_spline(xx, spli, kind='cubic') #function from spline interpolation
		aTgas = []
		aRate = []
		for i in range(imax):
			Tmin = log10(1e3) #as log(T/K)
			Tmax = log10(1e8) #as log(T/K)
			if(Tmin>Tmax): sys.exit("ERROR: Tmin>Tmax! dE:"+str(dE)+" Tgas:"+str(Tmax) )
			Tgas = 1e1**(float(i)/(imax-1)*(Tmax-Tmin)+Tmin)
			#Tgas = kte*dE*1.57888e5 #floating is 1/kboltzman_Ry
			kte = Tgas/dE/1.57888e5
			#prepare xt
			if(ty in [1, 4]):
				xt = 1e0-log(cc)/log(kte+cc)
			elif(ty in [2,3,5,6]):
				xt = kte/(kte+cc)
			else:
				sys.exit("ERROR: unknown type "+str(ty))

			sups = ff(xt) #read interpolated value
			#prepare ups
			if(ty==1):
				ups = sups * log(kte+exp(1e0))
			elif(ty==2):
				ups = sups
			elif(ty==3):
				ups = sups/(kte+1e0)			
			elif(ty==4):
				ups = sups*log(kte+cc)			
			elif(ty==5):
				ups = sups/kte			
			elif(ty==6):
				ups = 1e1**sups
			else:
				sys.exit("ERROR: unknown type "+str(ty))
			ups = max(ups, 0e0)
			Te = Tgas * 8.6173324e-5 #K->eV
			TRy = Tgas * 1.3806488e-16 * 4.58742e10 #K->erg->Ry
			wi = gmult[levLow] #mult level for starting level
			wj = gmult[levUp]
			#excitation
			#rate = 8.629e-6/sqrt(Te) * ups / wi * exp(-dE/TRy) #cm3/s 8.63e-6, 8.01140884e-8
			#de-excitation
			rate = 8.629e-6 * ups / sqrt(Tgas) / wj  #cm3/s
			aTgas.append(format_double("%e" % log10(Tgas)))
			aRate.append(format_double("%e" % max(log10(rate+1e-40),-40)))

		allrate[str(levLow)+"->"+str(levUp)+"_"+collname] = {"aTgas":aTgas, "aRate":aRate}

		sTgas = ", ".join(aTgas)
		sRate = ", ".join(aRate)
		krate = "1d1**fspline((/"+sTgas+"/), (/"+sRate+"/), logTgas)"
		fout.write("\n"+collname+", "+str(levUp-1)+", "+str(levLow-1)+", "+krate+"\n")

	return allrate

fout = open(foutname,"w")

#o_4.elvlc  o_4.psplups  o_4.splups  o_4.wgfa
for fff in flist:
	print "*****************"
	pre = fff.split("/")[-1]
	fout.write("\n#############\n")
	fout.write("#"+pre+"\n")
	apre = pre.split("_")
	nion = int(apre[1])-1
	if(nion==0):
		post = ""
	if(nion==1):
		post = "+"
	elif(nion>1):
		post = "+"*(nion)
	else:
		post = ""

	metal = apre[0].upper()+post
	metal_base = metal.replace("+","").replace("-","").lower()
	if(not(metal_base in usemetal)):
		print "skipping "+metal+" since not present in the array ",usemetal
		continue
	print metal
	fout.write("metal:"+metal+"\n")
	elv = fff+"/"+pre+".elvlc" #energy file
	wgf = fff+"/"+pre+".wgfa" #radiative data
	spl = fff+"/"+pre+".splups" #collisional (e-) data
	pspl = fff+"/"+pre+".psplups" #collisional (H+) data
	rrpfile = fff+"/"+pre+".rrparams"
	drpfile = fff+"/"+pre+".drparams"

	#****ELVLC****
	if(not(file_exists(elv))):
		"WARNING: "+elv+"doesn't exist, skipping..."
		continue
	if(not(file_exists(wgf))):
		"WARNING: "+wgf+"doesn't exist, skipping..."
		continue

	print "reading "+elv
	nparts = 10 #default size for ELVLC file (0=automatic)
	fh = open(elv,"rb") #open file
	gmult = dict()
	energy_data = dict()
	maxlev = 0 #maximum number of levels found
	#loop on file
	fout.write("#level n: energy (K), degeneracy g\n")
	Eold = 0e0 #previous energy found
	for row in fh:
		srow = row.strip()
		if(srow==""): continue #skip blank
		if(srow[0]=="%"): break #if comment breaks
		arow = [x for x in srow.split(" ") if is_number(x)] #keep only columns with numbers
		if(arow[0]=="-1"): break #if end data breaks
		if(nparts==0): nparts = len(arow) #store number of columns
		if(len(arow)<nparts):
			print arow
			sys.exit("ERROR: problem with the number of columns for "+elv) #check if number of cols is ok
		Enrg = float(arow[7]) #(observed), arow[9] (theoretical), unit: Ry
		if(Enrg==0e0): Enrg = float(arow[9]) #if observed not available use theoretical
		if(Enrg<Eold): Enrg = float(arow[9]) #if lower use theoretical
		level = int(arow[0])
		#check if energy data are available
		if(Enrg==0e0 and level>1):
			print srow
			sys.exit("ERROR: reading non strictly positive energy for a non-ground level!")
		#check energy
		#if(Enrg<Eold):
		#	print srow
		#	print "old","new",Eold,Enrg
		#	print "WARNING: Found energy lower than the previous level!"
			
		gmult[level] = int(arow[5]) #multeplicity
		Enrg_K = Enrg*1.57888e5
		if(level>levLimit): continue
		energy_data[level] = Enrg_K #K
		maxlev = max(maxlev,level)
		fout.write("level "+str(level-1)+": "+str("%e" % (Enrg*1.57888e5))+", "+str(gmult[level])+"\n")
		if(plotpre in pre): l = plt.axhline(y=level,ls="--")
		Eold = Enrg
	print "levels found:",maxlev
	print "max energy (K):",("%e" % energy_data[maxlev])

		

	#****WGFA****
	print "reading "+wgf
	fh = open(wgf,"rb") #open file
	fout.write("\n#up -> low, Aij (1/s), wavelength (angstrom)\n")
	itrans = 0
	for row in fh:
		srow = row.strip()
		if(srow==""): continue #skip blank
		if(srow[0]=="%"): break #if comment breaks
		if(srow=="-1"): break
		#arow = [x for x in srow.split(" ") if is_number(x)] #keep only columns with numbers
		arow = [row[j*5:(j+1)*5] for j in range(2)]
		arow += [x for x in row[10:].strip().split(" ") if is_number(x)] #keep only columns with numbers
		levLow = int(arow[0])
		levUp = int(arow[1])
		if(levLow>levLimit or levUp>levLimit): continue
		Aij = float(arow[4])
		wvl = float(arow[2])
		deltaE = 0e0
		if(abs(energy_data[levUp]-energy_data[levLow])>0e0):
			deltaE = 6.6260755e-27*2.99792458e10/abs(energy_data[levUp]-energy_data[levLow])/1.3806488e-16 * 1e8 #AA
		#print deltaE,wvl
		if(plotpre in pre): l = plt.axvline(x=itrans, ymin=float(levLow-1.)/(maxlev+1), ymax=float(levUp-1.)/(maxlev+1))
		fout.write(str(levUp-1)+" -> "+str(levLow-1)+", "+str("%e" % Aij)+", "+str("%e" % wvl)+"\n")
		itrans += 1
	print "transitions found:",itrans
		
	#read and prepare de-excitation rates for e- and H+ colliders
	arates = read_kex(spl,"e")
	read_kex(pspl,"H+")

	#print coefficient if needed
	#if(metal=="C++++"):
	#	for k,v in arates.iteritems():
	#		for ii in range(len(v["aTgas"])):
	#			print k,v["aTgas"][ii],v["aRate"][ii]
	#		print

	fout.write("\nendmetal\n")

plt.show()
print
print "*********************"
print "Stored in "+foutname
print "You can use it with KROME with the option -coolFile=path/"+foutname
print " where the path is relative to the ./krome script."


