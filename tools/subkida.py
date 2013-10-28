#convert from KIDA file to KROME network (including a set of option to chose a subset)

#input filename (provided by KIDA website)
fname = "kida_demo.dat" #KIDA filename
#output filename for KROME
foutname = "react_subkida" #output file

#several options to select your own reaction subset
avoid = [] #avoid atoms (can be empty)
use = ["e","C","H","N"] #use atoms (can be empty, e=electron)
maxatoms = 6 #maximum number of atoms
ions = True #use ions
anions = True #use anions
Tmin = -9999. #minimum temperature
Tmax = 99999. #maximum temperature
exclude = ["C8"] #species to exclude (can be empty)
#Recommendation is the recommendation given by experts in KIDA. (from KIDA help)
#0 means that the value is not recommended. 
#1 means that there is no recommendation (experts have not looked at the data). 
#2 means that the value has been validated by experts over the Trange
#3 means that it is the recommended value over Trange.
recom = [1,2,3] #recomandations to include (see above)

#Processes (selected using fitting formula)
#1: Cosmic-ray ionization (direct and undirect)
#2: Photo-dissociation (Draine)
#3: Kooij 
#4: ionpol1 
#5: ionpol2 
processes = [3,4,5] #processes included (see above)

#variables for cosmic rays and photochemistry
CRvar = "crflux" #name of the CR flux variable
Avvar = "Av" #name of the Av variable


#******************************************************************
#******************************************************************
#******************************************************************
# MODIFY UNDER THIS LINE ONLY IF YOU KNOW WHAT TO DO!
#******************************************************************
#******************************************************************
#******************************************************************

use = [x.upper() for x in use]
avoid = [x.upper() for x in avoid]
exclude = [x.upper() for x in exclude]
recom = [int(x) for x in recom]
processes = [int(x) for x in processes]

import sys

#function to check if argument is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#************************
#molecule class
class mol():
	name = "" #name
	atoms = [] #list of atoms
	natoms = 0 #count atoms
	charge = 0 #charge
	emol = [] #exploded molecule HCO2H = [H,H,C,O,O]
	#method to parse the name
	def parse(self):
		rms = ["c-","l-","(1)","+","-","(2D)","(1D)"]
		adic = [str(i) for i in range(30)] #numbers
		adic += ["E","H","HE","LI","C","N","O","F","NE","NA","MG","AL","SI","P","S","CL","CA","K","MN","FE","NI"] #elements
		adic = sorted(adic, key=lambda x:len(x),reverse=True) #sort by lenght
		name = self.name #copy name locally
		cname = name
		for x in rms:
			cname = cname.replace(x,"")
		cname = cname.upper() #copy name and remove charge
		amol = []
		atoms = []
		#loop to find multiple objects
		for j in range(10):
			#loop on atoms
			for a in adic:
				#serach for element
				if(a in cname):
					amol.append([cname.index(a),a]) #append position and name
					cname = cname.replace(a,"X"*len(a),1) #replace element with XX
					if(not(a in atoms) and not(is_number(a))): atoms.append(a) #store atom
		amol = sorted(amol,key=lambda x: x[0]) #sort by position
		emol = [] #exploded molecule
		if(cname.replace("X","").strip()!=""):
			print name,cname
			sys.exit()
		i = 0
		while(True):
			mult = 1 #multeplicity (e.g. H2 => mult=2)
			if(i>=len(amol)): break #end loop
			#if number skip
			if(is_number(amol[i][1])):
				i += 1
				continue
			#if next exists
			if(i<len(amol)-1):
				#if next is number => is multeplicity
				if(is_number(amol[i+1][1])): mult = int(amol[i+1][1])
			emol += ([amol[i][1]] * mult) #add to exploded
			i += 1 
		self.emol = emol #exploded
		self.atoms = atoms #atoms
		self.natoms = len(emol) #count atoms
		self.charge = name.count("+")-name.count("-") #charge

#**************************************
#**************************************
#**************************************

#Format : 3(a11) 1x 5(a11) 1x 3(e10.3 1x) 2(e8.2) 1x a4 i3 2(i7) i3 4x 2(i2) 1x i2 
#Reactants   Products  alpha  beta  gamma  F g Type_of_uncertainty   itype    Tmin  Tmax  Formula   Number     Number_of_(alpha, beta, gamma)
# Recommendation

print "************************"
print "*******SUB KIDA!********"
print "************************"

print "reading file "+fname+", wait..."

okcount = totcount = 0
nhist = dict()
idxdic = dict()
fh = open(fname,"rb")
fout = open(foutname,"w")
fmt = [11]*3 + [1] + 5*[11] +[1] + 3*[11] + [8,9] + [1,4,3] + 2*[7] + [3,6,3,2]
keys = ["R"+str(i) for i in range(3)] + ["x"] + ["P"+str(i) for i in range(5)] +["x"] + ["a","b","c"] + ["F","g"] + ["x","unc","type"] + ["tmin","tmax"] + ["formula","num","subnum","recom"]
reacts = []
fout.write("#reaction subset from KIDA created with KROME (see Grassi+2013)\n")
if(len(avoid)>0): fout.write("#avoid="+(",".join(avoid))+"\n")
if(len(use)>0): fout.write("#use="+(",".join(use))+"\n")
if(len(exclude)>0): fout.write("#exclude="+(",".join(exclude))+"\n")
fout.write("#recom="+(",".join([str(x) for x in recom]))+"\n")
fout.write("#processes="+(",".join([str(x) for x in processes]))+"\n")
if(maxatoms<100): fout.write("#maxatoms="+str(maxatoms)+"\n")
if(not(ions)): fout.write("#no ions\n")
if(not(anions)): fout.write("#no anions\n")
fout.write("#Tmin="+str(Tmin)+" Tmax="+str(Tmax)+"\n")
fout.write("@format:idx,R,R,R,P,P,P,P,P,Tmin,Tmax,rate\n")
rems = ["Photon","CRP","CRPHOT","CR",""]
for row in fh:
	totcount += 1
	srow = row.strip()
	if(srow==""): continue #skip empty lines
	if(srow[0]=="#"): continue #skip comments
	p = 0
	arow = dict()
	for i in range(len(fmt)):
		arow[keys[i]] = srow[p:p+fmt[i]].strip()
		p += fmt[i]
	RR = ",".join([arow["R"+str(i)] for i in range(3)])
	PP = ",".join([arow["P"+str(i)] for i in range(5)])
	TT = arow["tmin"]+","+arow["tmax"]

	#Formula is a number that referes to the formula needed to compute the rate coefficient of the reaction. 
	#see http://kida.obs.u-bordeaux1.fr/help
	#1: Cosmic-ray ionization (direct and undirect)
	#2: Photo-dissociation (Draine)
	#3: Kooij 
	#4: ionpol1 
	#5: ionpol2 
	#6: Troe fall-off (NOT SUPPORTED!)
	arow["formula"] = int(arow["formula"])
	if(arow["formula"]==1):
		KK = arow["a"]+"*"+CRvar
	elif(arow["formula"]==2):
		KK = arow["a"]
		if(float(arow["c"])!=0e0): KK += "*exp(-"+arow["c"]+"*"+Avvar+")"
	elif(arow["formula"]==3):
		KK = arow["a"]
		if(float(arow["b"])!=0e0): KK += "*(T32)**("+arow["b"]+")"
		if(float(arow["c"])!=0e0): KK += "*exp(-"+arow["c"]+"*invT)"
	elif(arow["formula"]==4):
		KK = arow["a"]
		if(float(arow["b"])!=1e0): KK += "*"+arow["b"]
		gpart = ""
		if(float(arow["c"])!=0e0): gpart = "+ 0.4767d0*"+arow["c"]+"*sqrt(3d2*invT)"
		KK += "*(0.62d0 "+gpart+")"
	elif(arow["formula"]==5):
		KK = arow["a"]
		if(float(arow["b"])!=1e0): KK += "*"+arow["b"]
		gpart = ""
		if(float(arow["c"])!=0e0): 
			gpart = "+ 0.0967d0*"+arow["c"]+"*sqrt(3d2*invT) + "
			gpart += arow["c"]+"**2*28.501d0*invT"
			KK += "*(1d0 "+gpart+")"
	elif(arow["formula"]==6):
		continue #WARNING: 3body not supported!
	else:
		print "ERROR: Formula not found!",arow["formula"]
	KK = KK.replace("--","+").replace("++","+").replace("-+","-").replace("+-","-")

	krow = (RR+","+PP+","+TT+","+KK+"\n")
	for x in rems:
		krow = krow.replace(x,"")
	
	ok = True
	#print arow
	if(not(arow["formula"] in processes)): ok = False
	if(not(int(arow["recom"]) in recom)): ok = False 
	if(float(arow["tmin"])<Tmin): ok = False
	if(float(arow["tmax"])>Tmax): ok = False
	if(not(ok)): continue
	for i in range(3):
		if(not(ok)): break
		mymol = mol()
		R = arow["R"+str(i)]
		if(R in rems): continue
		if(R in exclude): ok = False
		mymol.name = R
		mymol.parse()
		if(mymol.natoms>maxatoms): ok = False
		if("+" in R and not(ions)): ok = False
		if("-" in R and not(anions)): ok = False
		for x in avoid:
			if(x in mymol.atoms): ok = False
		if(len(use)>0):
			for x in mymol.atoms:
				if(not(x in use)): ok = False

	for i in range(5):
		if(not(ok)): break
		mymol = mol()
		P = arow["P"+str(i)]
		if(P in rems): continue
		if(P in exclude): ok = False
		mymol.name = P
		mymol.parse()
		if(mymol.natoms>maxatoms): ok = False
		if("+" in P and not(ions)): ok = False
		if("-" in P and not(anions)): ok = False
		for x in avoid:
			if(x in mymol.atoms): ok = False
		if(len(use)>0):
			for x in mymol.atoms:
				if(not(x in use)): ok = False
	if(not(ok)): continue
	if(arow["num"] in nhist): 
		nhist[arow["num"]].append(arow["subnum"])
		idxdic[arow["num"]] = okcount+1
	else: nhist[arow["num"]] = [arow["subnum"]]

	okcount += 1
	fout.write(str(okcount)+","+krow)
	reacts.append(arow)

fout.close()
fh.close()


#FINAL OUTPUT
print "Rections included:", okcount
if(okcount==0):
	print "********** WARNING! **********"
	print "No reactions match your criteria!!!"
	sys.exit()
print "Total reactions:", totcount
print "File written in:",foutname
multi = [str(k)+"("+str(idxdic[k])+")" for (k,v) in nhist.iteritems() if len(v)>1 ]
if(len(multi)>1):
	print "-------------------"
	print "The following reactions have multiple values:"
	print " "+(", ".join(multi))
	print " as KIDA_index (KROME_index)"
	print " please check!"
	print "-------------------"

print "Everything done! Bye!"
print


