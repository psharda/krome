#!/usr/bin/python
#!python

#THIS SCRIPT PREPARES COOLING TABLES WITH SELF-SHIELDING
# BY READING DATA FROM LAMDA-like FILES
#written by Tommaso Grassi and the KROME team (Apr3,2014)

from scipy.optimize import fsolve
from scipy import interpolate
import math,sys

#filename
fname = "coolCO.dat"

#amount of colliders
xH2 = 1e-2
xcoll = {"H2pa": 0.25e0*xH2, "H2or": 0.75*xH2}

#temperature range
tmin = math.log10(2e1)
tmax = math.log10(2900e0)

#temperature steps
imax = 30

#velocity gradient cm/s/cm
absdvdz = abs(1e-11)

############################################
# PLEASE DO NOT MODIFY UNDER THIS LINE IF
# YOU ARE NOT SURE OF WHAT YOU ARE DOING!
############################################

#some constants
kboltzmann = 1.380658e-16 #erg/K
hplanck = 6.6260755e-27 #erg*s
clight = 2.99792458e10 #cm/s

#define some dictionaries and lists
energy = dict() #energy in K, key=level number
g = dict() #multeplicity, key=level number
nu = dict() #frequency in 1/s, key=level number
A = dict() #Aij coefficient, key=(i,j) tupla
B = dict() #B coefficient, key=(i,j) tupla
Cfit = dict() #Collisional coeff, key=collider; this is a dicionary of dictionaries
colliders = [] #list of colliders

print "**************************"
print "        TABCOOL!"
print "**************************"
print "Reading from "+fname

#read data from file
fh = open(fname)

#loop on file
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	if(srow[0]=="#"): continue
	srow = srow.split("#")[0]
	srow = srow.split("//")[0]

	#read and store level data: number, energy, multeplicity
	#compute also frequency
	if("level" in srow): 
		arow = srow.split(":")
		level_number = int(arow[0].replace("level","").strip()) #level number (zero-based)
		(level_energy, level_g) = [x.strip() for x in arow[1].split(",")]
		energy[level_number] = float(level_energy) #energy in K
		nu[level_number] = float(level_energy) * kboltzmann / hplanck #frequency in 1/s
		g[level_number] = float(level_g) #level multeplicity

	#read and store spontaneus transitions: upper level, lower level, Aij coefficient
	if(len(srow.split(","))==3):
		(level_up, level_low, level_A) = [x.strip() for x in srow.split(",")]
		level_up = int(level_up) #upper level
		level_low = int(level_low) #lower level
		level_A = float(level_A) #Aij coefficient 1/s
		A[(level_up, level_low)] = level_A #store Aij
		#compute Bji and Bij from Aij
		B[(level_up, level_low)] = .5e0 * level_A * clight**2 / hplanck / (nu[level_up] - nu[level_low])**3
		B[(level_low, level_up)] = B[(level_up, level_low)] * g[level_up] / g[level_low]
	
	#retrive collisional data from flin F90-style function
	# data are: collider name, upper level, lower level, rate in F90-style
	if("flin" in srow):
		(collider, collider_up, collider_down, rate) = srow.split(",",3)[:4]
		#explode flin to obtain Trange and corresponding function values
		arate = rate.replace("d","e").replace("(","").replace(")","").replace(", Tgas","").replace("flin","").split("/, /")
		arate = [x.replace("/","").strip() for x in arate]
		Trange = [float(x) for x in arate[0].split(",")] #store Trange values in a list
		vrate = [float(x) for x in arate[1].split(",")] #store coefficient values in a list
		if(len(Trange)!=len(vrate)): sys.exit("ERROR: Trange must be equal to vrate!")
		#store collider names and define the dictionary of Cij for the given collider
		if(not(collider in colliders)):
			colliders.append(collider)
			Cfit[collider] = dict()
		#store Trange and vrate in the dictionary
		Cfit[collider][(int(collider_up), int(collider_down))] = {"T":Trange, "val":vrate}

#function to compute tau according to NK93
def ftau(lBji,lBij,dv,xj,xi):
	return hplanck*clight*0.25/math.pi/dv*(xj*lBji-xi*lBij)

#function to beta according to NK93
def fbeta(lBji,lBij,dv,xj,xi):
	mytau = ftau(lBji,lBij,dv,xj,xi)
	return 1e0/(1e0+3e0*mytau)

#define system of equations as in NK93
def eqs(x,args):
	#read data from additional arguments
	dvdz = args["dvdz"]
	lA = args["A"]
	lB = args["B"]
	lC = args["C"]
	lcolliders = args["colliders"]
	lxcoll = args["xcoll"]
	nlev = args["nlev"]
	eq = [0e0 for i in range(nlev)] #initialize equations
	#loop on levels to prepare equations
	for i in range(nlev):
		#prepare first term (j->i)
		p1 = 0e0
		#loop on levels
		for j in range(nlev):
			if(i==j): continue #i=j is not a transition!
			#use only available transitions
			try:
				AA1 = lA[(j,i)]
				AA1 *= fbeta(lB[(i,j)],lB[(j,i)],dvdz,x[i],x[j]) #add shielding
			except:
				AA1 = 0e0 #otherwise transition is zero
			CC = 0e0
			#sum up (de)excitations from colliders
			for coll in lcolliders:
				CC += lC[coll][(j,i)] * lxcoll[coll]
			#put everything together
			p1 += x[j] * (AA1 + CC)

		#prepare second term (i->j), see comments above
		p2 = 0e0
		for j in range(nlev):
			if(i==j): continue #i=j is not a transition!
			try:
				AA2 = lA[(i,j)]
				AA2 *= fbeta(lB[(j,i)],lB[(i,j)],dvdz,x[j],x[i])
			except:
				AA2 = 0e0

			CC = 0e0
			for coll in lcolliders:
				CC += lC[coll][(i,j)] * lxcoll[coll]

			p2 += (AA2 + CC)

		#build equation for the ith level
		eq[i] = p1 - x[i] * p2
	
	#first equation is replaced by continuity
	eq[0] = sum(x) - 1e0
	return eq

#loop on temperatures
print "starting calculations..."
print "Tgas, cooling (erg/s/cm3)"
for ii in range(imax):
	Tgas = 1e1**(ii*(tmax-tmin)/imax+tmin) #compute Tgas
	C = dict() #init collision dictionary
	#loop on colliders
	for coll in colliders:
		C[coll] = dict()
		#loop on transitions for the given collider
		for k,v in Cfit[coll].iteritems():
			#interpolate Cij for Tgas
			f = interpolate.interp1d(Cfit[coll][k]["T"], Cfit[coll][k]["val"])
			C[coll][k] = f(Tgas) #store Cij
			deltaE = energy[k[0]]-energy[k[1]] #already K
			#compute and store Cji from Cij
			C[coll][(k[1],k[0])] = C[coll][k] * g[k[0]]/g[k[1]] * math.exp(-deltaE/Tgas)


	#prepare additional arguments
	myargs = {"A":A, "B":B, "C":C,"nlev":len(energy),"dvdz":absdvdz,"colliders":colliders,"xcoll":xcoll}

	#initial guess
	x0 = [1e0/len(energy) for i in range(len(energy))]
	
	#solve the (non)linear system
	res = fsolve(eqs, x0,args=myargs,full_output=True)
	if(res[2]!=1): sys.exit("ERROR: "+res[3]) #check for error
	xx = res[0] #copy solution
	#compute cooling as x(j)*Aji*betaji (erg/s/cm3)
	cool = 0e0
	for k,v in A.iteritems():
		kup = k[0]
		klow = k[1]
		beta = fbeta(B[(klow,kup)],B[(kup,klow)],myargs["dvdz"],xx[klow],xx[kup])
		cool += A[k] * xx[kup] * (energy[kup] - energy[klow]) * kboltzmann
	print Tgas,cool

print "DONE!"




