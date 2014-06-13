#!/usr/bin/python
#!python

import sys
from os.path import exists as file_exists
from math import log,log10,exp,sqrt
from scipy.interpolate import interp1d as interp_spline
from os import listdir
from os.path import isfile, join,isdir

foutname = "coolChianti.dat"
path = "chianti/"
flist = []
dirs = [ f for f in listdir(path) if isdir(join(path,f)) ]
for dr in dirs:
	mypath = path+dr+"/"
	dirs2 = [f for f in listdir(mypath) if isdir(join(mypath,f)) ]
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


fout = open(foutname,"w")

#o_4.elvlc  o_4.psplups  o_4.splups  o_4.wgfa
for fff in flist:
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
	fout.write("metal:"+metal+"\n")
	elv = fff+"/"+pre+".elvlc" #energy file
	wgf = fff+"/"+pre+".wgfa" #radiative data
	spl = fff+"/"+pre+".splups" #collisional (e-) data
	pspl = fff+"/"+pre+".psplups" #collisional (H+) data

	#****ELVLC****
	if(not(file_exists(elv))): sys.exit("ERROR: "+elv+" does not exists!")
	print "reading "+elv
	nparts = 10 #default size for ELVLC file (0=automatic)
	fh = open(elv,"rb") #open file
	gmult = dict()
	#loop on file
	fout.write("#level n: energy (K), degeneracy g\n")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue #skip blank
		if(srow[0]=="%"): break #if comment breaks
		arow = [x for x in srow.split(" ") if is_number(x)] #keep only columns with numbers
		if(arow[0]=="-1"): break #if end data breaks
		if(nparts==0): nparts = len(arow) #store number of columns
		if(len(arow)!=nparts): sys.exit("ERROR: problem with the number of columns for "+elv) #check if number of cols is ok
		Enrg = float(arow[7]) #(observed), arow[8] (theoretical), unit: Ry
		level = int(arow[0])
		gmult[level] = int(arow[5]) #multeplicity
		fout.write("level "+str(level-1)+": "+str("%e" % (Enrg*1.57888e5))+", "+str(gmult[level])+"\n")

		

	#****WGFA****
	if(not(file_exists(wgf))): sys.exit("ERROR: "+wgf+" does not exists!")
	print "reading "+wgf
	fh = open(wgf,"rb") #open file
	fout.write("\n#Aij (1/s)\n")
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
		Aij = float(arow[4])
		fout.write(str(levUp-1)+", "+str(levLow-1)+", "+str("%e" % Aij)+"\n")

		

	#****SPLUPS****
	#method: Burgess+Tully 1992 A&A
	if(not(file_exists(spl))): sys.exit("ERROR: "+spl+" does not exists!")
	print "reading "+spl
	fh = open(spl,"rb") #open file
	fout.write("\n#collider, level_up, level_down, rate\n")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue #skip blank
		if(srow[0]=="%"): break #if comment breaks
		if(srow=="-1"): break #if end data breaks
		arow = [row[j*3:(j+1)*3] for j in range(5)]
		for j in range(len(row[15:])/10):
			block = [row[15+j*10:15+(j+1)*10].strip()]
			if(block!=[""]): arow += block
		#arow += [x for x in row[15:].strip().split(" ") if is_number(x)] #keep only columns with numbers
		levLow = int(arow[2])
		levUp = int(arow[3])
		ty = int(arow[4]) #type
		dE = float(arow[6]) #Ry
		cc = float(arow[7]) #scaling factor
		spli = [float(x) for x in arow[8:]] #read spline
		xx = [float(x)/(len(spli)-1) for x in range(len(spli))] #normalized x value for spline
		imax = int(1e1) #number of Tgas bins
		ff = interp_spline(xx, spli, kind='cubic') #function from spline interpolation
		aTgas = []
		aRate = []
		for i in range(imax):
			Tmin = 3e0 #as log(T/K)
			Tmax = log10(400000) #as log(T/K)
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
			Te = Tgas * 8.6173324e-5 #K->eV
			TRy = Tgas * 1.3806488e-16 * 4.58742e10 #K->erg->Ry
			wi = gmult[levLow] #mult level for starting level
			rate = 8.629e-6/sqrt(Te) * ups / wi * exp(-dE/TRy) #cm3/s 8.63e-6, 8.01140884e-8
			aTgas.append(format_double("%e" % Tgas))
			aRate.append(format_double("%e" % max(rate,1e-40)))

		sTgas = ", ".join(aTgas)
		sRate = ", ".join(aRate)
		krate = "flin((/"+sTgas+"/), (/"+sRate+"/), Tgas)"
		fout.write("\ne, "+str(levUp-1)+", "+str(levLow-1)+", "+krate+"\n")
	fout.write("\nendmetal\n")


print "Stored in "+foutname
print "You can use it with KROME with the option -coolFile=path/"+foutname
print " where the path is relative to the ./krome script."


