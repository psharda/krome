fh = open("photo.dat","rb")
fout = open("kromeauto.dat","w")

def format_double(snum):
	snum = str(snum).lower()
	#format string number to F90 double
	if("d" in snum): return snum
	if("e" in snum): return snum.replace("e","d")	
	return snum+"d0"

zdic = {"E":0,
	"H":1,
	"HE":2,
	"LI":3,
	"BE":4,
	"B":5,
	"C":6,
	"N":7,
	"O":8,
	"F":9,
	"NE":10,
	"NA":11,
	"MG":12,
	"AL":13,
	"SI":14,
	"P":15,
	"S":16,
	"CL":17,
	"AR":18,
	"K":19,
	"CA":20,
	"SC":21,
	"TI":22,
	"V":23,
	"CR":24,
	"MN":25,
	"FE":26,
	"CO":27,
	"NI":28,
	"CU":29,
	"ZN":30
}

zz = dict()
for k,v in zdic.iteritems():
	zz[v] = k

for row in fh:
	srow = row.strip()	
	arow = [x for x  in srow.split(" ") if x!=""]
	el = zz[int(arow[0])].title()
	elo = el
	dz = int(arow[0]) - int(arow[1])
	eli = el+"+"
	if(dz==1):
		el += "+"
		eli = el+"2" 
	elif(dz>1):
		el += "+" +str(dz)
		eli = elo+"+"+str(dz+1)
	#(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
	Eth = arow[2]
	Emax = arow[3]
	fout.write("#photoionization rate from verner+96\n")
	fout.write("@type: photoion\n")
	fout.write("@reacts: "+el+"\n")
	fout.write("@prods: "+eli+", E\n")
	fout.write("@limits: "+format_double(Eth)+", "+format_double(Emax)+"\n")
	xx = [format_double(x) for x in arow[4:]]
	fout.write("@rate: sigma_v96(energy_ev, "+(", ".join(xx))+")\n\n")

fh.close()
fout.write("\n######################################\n")
fout.write("######################################\n")


fh = open("table1.dat","rb")
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	arow = [x for x  in srow.split(" ") if x!=""]

	el = zz[int(arow[1])].title()
	elo = el
	eli = el+"+"
	dz = int(arow[1]) - int(arow[2])
	if(dz==1):
		el += "+"
		eli = el+"2" 
	elif(dz>1):
		el += "+" +str(dz)
		eli = elo+"+"+str(dz+1)

	fout.write("#radiative recombination rate from verner+96\n")
	fout.write("@type: radrec\n")
	fout.write("@reacts: "+eli+", E\n")
	fout.write("@prods: "+el+"\n")
	fout.write("@limits: 3d0, 1d10\n")
	xx = [format_double(x) for x in arow[3:]]
	fout.write("@rate: radrec_v96(Tgas, "+(", ".join(xx))+")\n\n")


fh.close()
fout.write("\n######################################\n")
fout.write("######################################\n")


fh = open("cfit.dat","rb")	
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	srow = srow[6:]
	arow = [x for x  in srow.split(" ") if x!=""]
	print arow
	el = zz[int(arow[0])].title()
	elo = el
	eli = el+"+"
	dz = int(arow[0]) - int(arow[1])
	if(dz==1):
		el += "+"
		eli = el+"2" 
	elif(dz>1):
		el += "+" +str(dz)
		eli = elo+"+"+str(dz+1)

	print el,eli,dz
	fout.write("#collisional ionization rate from verner+96\n")
	fout.write("@type: colion\n")
	fout.write("@reacts: "+el+", E\n")
	fout.write("@prods: "+eli+",E, E\n")

	Emin = 1.1604519e4*float(arow[7]) #eV->K
	Emax = 1.1604519e4*float(arow[8]) * 1e3 #keV->K
	fout.write("@limits: "+format_double('%e' % Emin)+", "+format_double('%e' % Emax)+"\n")
	xx = [format_double(x) for x in arow[2:7]]
	fout.write("@rate: colion_v96(Tgas, "+(", ".join(xx))+")\n\n")




