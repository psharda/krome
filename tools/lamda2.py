#!/usr/bin/python
#!python

#THIS PYTHON SCRIPT IS A PART OF KROME AND ALLOWS TO CONVERT THE 
# SPECTROSCOPIC DATA FROM THE LAMDA DATABASE TO THE KROME FORMAT
# FOR METAL/MOLECULAR COOLING
# WRITTEN by Tommaso Grassi and the KROME group.
import sys,os

########################################
# PLEASE DO NOT MODIFY IF YOU ARE NOT 
# SURE OF WHAT YOU ARE DOING!
########################################

#name = "co.dat" #input file
#outfile = "out.dat" #output file

#check command-line arguments and store them
if(len(sys.argv)<3):
    sys.exit("Usage: %s INPUT OUTPUT" % sys.argv[0])

if(not os.path.isfile(sys.argv[1].strip())):
    sys.exit("ERROR: input file %s was not found!" % sys.argv[1])

fname = sys.argv[1]
outfile = sys.argv[2]
if(fname==outfile):
    sys.exit("ERROR: input and output are the same file!")

#function to read a line, strip and pack into a list
def read_line(line):
	line = line.replace("\t"," ")
	while("  " in line):
		line = line.replace("  "," ")
	aarow = line.strip().split(" ")
	return [ael.strip() for ael in aarow]

#star reading
# simple script, just reads lines knowing their position within the file
print "Now reading from "+fname
fh = open(fname,"rb")
data = dict()
rows = [x.strip() for x in fh]
fh.close()

#MOLECULE
rdx = 1
data["molecule"] = read_line(rows[rdx])[0]

#WEIGHT
rdx += 2
data["weight"] = float(rows[rdx])

#NLEVELS
rdx += 2
data["nlevels"] = int(rows[rdx])

#LEVELS
data["levels"] = []
rdx += 1
for i in range(data["nlevels"]):
	rdx += 1
	arow = read_line(rows[rdx])
	data["levels"].append({"level":arow[0], "energy":arow[1], "weight":arow[2], "J":arow[3]})

#NUMBER OF RADIATIVE TRANSITIONS
rdx += 2
data["nrad"] = int(rows[rdx])

#READ RADIATIVE TRANSITIONS
data["rads"] = []
rdx += 1
for i in range(data["nrad"]):
	rdx += 1
	arow = read_line(rows[rdx])
	tdic = {"irad":int(arow[0]), "up":arow[1], "low":arow[2], "Aij":arow[3], "nu":float(arow[4]), "energy":float(arow[5])}
	data["rads"].append(tdic)

#COLLISIONAL PARTNERS
rdx += 2
data["npartners"] = int(rows[rdx])

#READ PARTNERS
data["colls"] = dict()
data["partners"] = []
isOrthoPara = False
for j in range(data["npartners"]):
	#READ COLLIDER NAME AND INFO
	rdx += 2
	arow = read_line(rows[rdx])
	#according to LAMDA: collision partner ID and reference.
	# Valid identifications are: 1=H2, 2=para-H2, 3=ortho-H2, 4=electrons, 5=H, 6=He.
	valid_id = ["","H2","H2pa","H2or","e","H","He"]
	partner_id = int(arow[0])
	if(partner_id>len(valid_id)-1 or partner_id==0): sys.exit("ERROR: unknown partner id "+str(partner_id))
	#partner_name = arow[1].replace(data["molecule"]+"-","").replace(":","")
	partner_name = valid_id[partner_id]
	if(partner_id==2 or partner_id==3):
		isOrthoPara = True
	data["colls"][partner_name] = dict()
	data["colls"][partner_name]["info"] = rows[rdx]
	data["colls"][partner_name]["isOrthoPara"] = isOrthoPara
	data["partners"].append(partner_name)
	#NUMBER OF COLLISIONAL TRANSITIONS
	rdx += 2
	data["colls"][partner_name]["ncolls"] = int(rows[rdx]) 
	#NUMBER OF temperatures
	rdx += 2
	data["colls"][partner_name]["ntemps"] = int(rows[rdx]) 
	#TEMPERATURES
	rdx += 2
	data["colls"][partner_name]["temps"] = read_line(rows[rdx]) 
	#READ TRANSITIONS ! TRANS + UP + LOW + RATE COEFFS(cm^3 s^-1)
	data["colls"][partner_name]["colls"] = []
	rdx += 1
	for i in range(data["colls"][partner_name]["ncolls"]):
		rdx += 1
		arow = read_line(rows[rdx])
		tdic = {"icoll":int(arow[0]), "up":arow[1], "low":arow[2], "rates":arow[3:]}
		data["colls"][partner_name]["colls"].append(tdic)

#PRINT DATA IN KROME FORMAT
fout = open(outfile,"w")
fout.write("############\n")
fout.write("metal:"+data["molecule"]+"\n")
fout.write("#level n: energy (K), degeneracy g\n")
for i in range(len(data["levels"])):
	nlev = str(int(data["levels"][i]["level"])-1)
	elev = str(float(data["levels"][i]["energy"])*1.42879)
	glev = data["levels"][i]["weight"]
	fout.write("level "+nlev+": "+elev+", "+glev+"\n")
fout.write("\n")
fout.write("#Aij (1/s)\n")
for i in range(len(data["rads"])):
	aup = int(data["rads"][i]["up"])-1
	alow = int(data["rads"][i]["low"])-1
	fout.write(str(aup)+" -> "+str(alow)+", "+data["rads"][i]["Aij"]+"\n")

#function to format numbers to F90 format (x.xdx)
def fmt_f90(xarg):
	xarg = xarg.lower()
	if("e" in xarg): 
		return xarg.replace("e","d")
	elif("d" in xarg): 
		return xarg
	else:
		return xarg + "d0"

#write the partners to the krome file
#note here that the function is employed to pass from lamda
# linear interpolation to krome linear interpolation.
for p in data["partners"]:
	fout.write("\n")
	fout.write("#collider, level_up, level_down, rate\n")
	pdata = data["colls"][p]
	if(pdata["isOrthoPara"]): fout.write("ortho/para: 3/1\n")
	for i in range(len(pdata["colls"])):
		myc = pdata["colls"][i]
		cup = str(int(myc["up"])-1)
		clow = str(int(myc["low"])-1)
		xvals = "(/"+(", ".join([fmt_f90(x)  for x in pdata["temps"]]))+"/)"
		yvals = "(/"+(", ".join([fmt_f90(x) for x in myc["rates"]]))+"/)"
		fout.write(p+", "+cup+", "+clow+", fspline("+xvals+", "+yvals+", Tgas)\n")

fout.write("end metal\n")

fout.close()

#RECURSIVE FUNCTION TO SHOW THE DATA STRUCTURE
def pd(d,depth):
	for (k,v) in d.iteritems():
		if(isinstance(v,dict) or isinstance(v,list)):
			ll = "x"+str(len(v))
		else: 
			ll = "{"+str(v)+"}"
		print ("---"*depth)+"> " + k,ll
		if(isinstance(v,dict)): pd(v,depth+1)

#print a pseudo-data structure, just to check that everything went smooth
print "Data structure:"
pd(data,0)
print "output written in "+outfile
print "done!"
