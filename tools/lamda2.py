#THIS PYTHON SCRIPT IS A PART OF KROME AND ALLOWS TO COVERT THE 
# SPCTROSCOPIC DATA FROM THE LAMDA DATABSE TO THE KROME FORMAT
# FOR METAL/MOLECULAR COOLING

fname = "hnc.dat" #input file
outfile = "out.dat" #output file


def read_line(line):
	line = line.replace("\t"," ")
	while("  " in line):
		line = line.replace("  "," ")
	aarow = line.strip().split(" ")
	return [ael.strip() for ael in aarow]


print "Now reading from "+fname
fh = open(fname,"rb")
data = dict()
rows = [x.strip() for x in fh]
fh.close()

#MOLECULE
rdx = 1
data["molecule"] = rows[rdx]

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

#NUMBER OF RADIATIVE TRANSTIONS
rdx += 2
data["nrad"] = int(rows[rdx])

#READ RADIATIVE TRANSITONS
data["rads"] = []
rdx += 1
for i in range(data["nrad"]):
	rdx += 1
	arow = read_line(rows[rdx])
	tdic = {"irad":int(arow[0]), "up":arow[1], "low":arow[2], "Aij":arow[3], "nu":float(arow[4]), "energy":float(arow[5])}
	data["rads"].append(tdic)

#COLLISIONAL PARTENERS
rdx += 2
data["npartners"] = int(rows[rdx])

#READ PARTNERS
data["colls"] = dict()
data["partners"] = []
for j in range(data["npartners"]):
	#READ COLLIDER NAME AND INFO
	rdx += 2
	arow = read_line(rows[rdx])	
	partner_name = arow[1].replace(data["molecule"]+"-","").replace(":","")
	data["colls"][partner_name] = dict()
	data["colls"][partner_name]["info"] = rows[rdx]
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
	fout.write(str(aup)+", "+str(alow)+", "+data["rads"][i]["Aij"]+"\n")

def fmt_f90(xarg):
	xarg = xarg.lower()
	if("e" in xarg): 
		return xarg.replace("e","d")
	elif("d" in xarg): 
		return xarg
	else:
		return xarg + "d0"

for p in data["partners"]:
	fout.write("\n")
	fout.write("#collider, level_up, level_down, rate\n")
	pdata = data["colls"][p]
	for i in range(len(pdata["colls"])):
		myc = pdata["colls"][i]
		cup = str(int(myc["up"])-1)
		clow = str(int(myc["low"])-1)
		xvals = "(/"+(", ".join([fmt_f90(x)  for x in pdata["temps"]]))+"/)"
		yvals = "(/"+(", ".join([fmt_f90(x) for x in myc["rates"]]))+"/)"
		fout.write(p+", "+cup+", "+clow+", flin("+xvals+", "+yvals+", Tgas)\n")
		#xfit = np.asarray([float(x) for x in pdata["temps"]])
		#yfit = np.asarray([float(x) for x in myc["rates"]])
		#p0 = sy.array([float(myc["rates"][0]), 1, .0001])
		#popt, pcov = curve_fit(ffit, xfit, yfit, p0)
		#print "****"
		#print popt,pcov
		#print i
		#for j in range(len(pdata["temps"])):
		#	xx = float(pdata["temps"][j])
		#	fgraph.write(str(i)+" "+str(pdata["temps"][j])+" "+myc["rates"][j]+" "+str(ffit(xx,popt[0],popt[1],popt[2]))+"\n")
fout.write("end metal\n")

fout.close()

def pd(d,depth):
	for (k,v) in d.iteritems():
		if(isinstance(v,dict) or isinstance(v,list)):
			ll = "x"+str(len(v))
		else: 
			ll = "{"+str(v)+"}"
		print ("---"*depth)+"> " + k,ll
		if(isinstance(v,dict)): pd(v,depth+1)

print "Data structure:"
pd(data,0)
print "output written in "+outfile
print "done!"
