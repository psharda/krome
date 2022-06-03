# KROME is a nice and friendly chemistry package for a wide range of
# astrophysical simulations. Given a chemical network (in CSV format)
# it automatically generates all the routines needed to solve the kinetic
# of the system, modelled as system of coupled Ordinary Differential
# Equations.
# It provides different options which make it unique and very flexible.
# Any suggestions and comments are welcomed. KROME is an open-source
# package, GNU-licensed, and any improvements provided by
# the users is well accepted. See disclaimer below and GNU License
# in gpl-3.0.txt.
#
# more details in http://kromepackage.org/
# also see https://bitbucket.org/krome/krome_stable
#
# more details in http://kromepackage.org/
# also see https://bitbucket.org/krome/krome_stable
#
# Written and developed by Tommaso Grassi
# tgrassi@nbi.dk,
# Starplan Center, Copenhagen.
# Niels Bohr Institute, Copenhagen.
#
# and Stefano Bovino
# stefano.bovino@uni-hamburg.de
# Hamburger Sternwarte, Hamburg.
#
# Contributors: J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
#  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
#  E.Tognelli
#
# KROME is provided "as it is", without any warranty.
# The Authors assume no liability for any damages of any kind
# (direct or indirect damages, contractual or non-contractual
# damages, pecuniary or non-pecuniary damages), directly or
# indirectly derived or arising from the correct or incorrect
# usage of KROME, in any possible environment, or arising from
# the impossibility to use, fully or partially, the software,
# or any bug or malefunction.
# Such exclusion of liability expressly includes any damages
# including the loss of data of any kind (including personal data)

# THIS FILE CONTAINS FUNCTIONS AND CLASSES FOR THE MAIN
# KROME PYTHON SCRIPT


import sys
import os
##################################
class molec:
	name = "" #molecule name
	charge = 0 #charge (0=neutral)
	mass = 0e0 #mass in g
	radius = 0 # 'radius' of species in cm
	zatom = 0 #atomic number (also for molecules)
	neutrons = 0 #total number of neutrons
	ename = "" #exploded name (e.g. H2C4=CCCCHH)
	fname = "" #f90 name (e.g. H+=Hj, D-=Dk)
	phname = "" #compact name for photoionizations (e.g. Fe++++=Fej4)
	coolname = "" #name for cooling, e.g. CIV
	fidx = "idx_" #index for fortran (e.g. H+=idx_Hj, D-=idx_Dk)
	is_atom = False #flag to identify atoms
	is_surface = False #flag for species on surface
	is_chemisorbed = False #flag for chemisorbed species
	is_clustarable = False #flag for species that can grow into larger clusters
	has_GFE_table = False #flag for having Gibbs free energy data
	# NOTE: Using JANAF tables is discouraged because it is unclear on the exact
	# calculations/data and it is incompatible with other databases like NASA
	# polynomial ones.
	# hasJanafThermoTable = False #flag for having thermo data from JANAF database
	Ebind_ice = Ebind_bare = 0e0 #binding energy on surface for ice and bare grain
	parentDustBin = 0 #for surface species: belongs to this dust bin (1-based)
	chempot = 0. #chemical potential (J/mol)
	poly1_nasa = [0.e0]*7 #nasa polynomials (usually 200-1000K)
	poly2_nasa = [0.e0]*7 #nasa polynomials (usually 1000-5000K)
	Tpoly_nasa = [0.e0]*3 #temperature limits
	poly1_nist = [0.e0]*7 #nist polynomials
	poly2_nist = [0.e0]*7 #nist polynomials
	Tpoly_nist = [0.e0]*3 #temperature limits
	idx = 0 #species index
	enthalpy = 0.e0 #enthalpy of formation
	atomcount = dict() #dictionary containin the count of atoms including zero (e.g H2O is {"H":2, "O":1, "C":0, ...})
	atomcount2 = dict() #dictionary containin the count of atoms without zero species (e.g H2O is {"H":2, "O":1})
	natoms = 0 #the number of atoms (e.g. diatomic=2)
	ve_vib = "__NONE__" #vibrational constant in K
	be_rot = "__NONE__" #rotational constant in K
	links = 0 #number of reactions involved
	nameLatex = "" #name in LaTeX format

	def __init__(self):
		self.poly1_nasa = [0.e0]*7
		self.poly2_nasa = [0.e0]*7
		self.Tpoly_nasa = [0.e0]*3
		self.poly1_nist = [0.e0]*7
		self.poly2_nist = [0.e0]*7
		self.Tpoly_nist = [0.e0]*3
		self.atomcount = dict()


##################################
class reaction:
	verbatim = "" #reaction written as string (e.g. A+B->C+D)
	reactants = [] #list of reactants (molec objects)
	products = [] #list of products (molec objects)
	Tmin = "0.d0" #min value of temperature range
	Tmax = "0.d0" #max value of temperature range
	TminOp = "" #min operator
	TmaxOp = "" #max operator
	pseudo_hash = "" #pseudo hash
	krate = "" #reaction rate written in F90 style
	idx = 0 #reaction index
	idxph = 0 #reaction index photorate (start from 1, 0=not a photorate)
	RHS = "" #ODE RHS in F90 style k(1)*n(2)*n(3)
	RHSvar = "" #ODE RHS as variable k1n2n3
	kphrate = None #photochemical rate
	dH = None #enthalpy of formation
	qeff = 0.e0 #effective Q value (for nuclear reactions)
	curlyR = [] #reactants curlyness
	curlyP = [] #products curlyness
	nuclearMult = "" #nuclear multeplicty factor 1/(n!)
	hasTlimitMax = hasTlimitMin = True #flag to determine the presence of Temperature limits
	group = "__DEFAULT__"
	canUseTabs = True #flag if this reaction can use tabs
	isStoreOnce = False #reaction rate is stored once and never evaluated again during fex call
	ifrate = "" #if condition on rate, e.g. if(Tgas>1d2):
	isCR = False #flag this reaction as CR
	isXRay = False #flag this reaction as XRay
	isSurface = False #flag this reaction as Surface reaction
	hasXsecFile = False #photo cross section is loaded from file
	xsecFile = "" #cross section file
	isAutoRev = False #this reaction is a reverse reaction in the automatic database
	photoForm = []
	photoDestroy = []

	#method: constructor to initialize lists
	def __init__(self):
		self.reactants = []
		self.products = []

	#method: build verbatim from reactants and products
	def build_verbatim(self):
		myr = []
		myp = []
		for r in self.reactants:
			if r.name != "dummy":
				myr.append(r.name)
		for p in self.products:
			if p.name != "dummy":
				myp.append(p.name)
		self.verbatim = " + ".join(myr)+" -> "+" + ".join(myp)

	#method: build photochemical rate
	def build_phrate(self, photoBlock=False):
		if "krome_kph_auto" not in self.krate and not photoBlock: return
		myr = self.reactants
		self.kphrate = self.krate.replace("krome_kph_auto=","")
		if "@xsecFile=" in self.krate:
			sidx = str(self.idx)
			args = "xsec"+sidx+"_val(:), xsec"+sidx+"_Emin,"
			#args += "xsec"+sidx+"_n,"
			args += "xsec"+sidx+"_idE, dshift("+self.reactants[0].fidx+")"
			self.kphrate = "xsec_interp(energyL, energyR, "+args+")"
			#self.xsecFile = self.krate.replace("@xsecFile=","").strip()
			if "SWRI" in self.krate.upper():
				RR = self.reactants[0].name
				PP = ("_".join([x.name for x in self.products]))
				self.xsecFile = "swri_"+RR+"__"+PP+".dat"
			if "LEIDEN" in self.krate.upper():
				RR = self.reactants[0].name
				PP = ("_".join([x.name for x in self.products]))
				self.xsecFile = "leiden_"+RR+"__"+PP+".dat"

			#replace with rate
			self.krate = self.krate.replace("@xsecFile=SWRI", "#XSECS#")
			self.krate = self.krate.replace("@xsecFile=LEIDEN", "#XSECS#")
			self.krate = self.krate.replace("#XSECS#", "photoBinRates("+str(self.idxph)+")")
		else:
			self.krate = "photoBinRates("+str(self.idxph)+")"

		#if(self.krate.strip()=="krome_kph_auto"):
		#	self.krate = "krome_kph_"+myr[0].phname
		#else:
		#	self.krate = "krome_kph_" + myr[0].phname.capitalize() + "R" + str(self.idx)

	#method: build RHS
	def build_RHS(self, useNuclearMult=False):
		if self.idx<=0:
			print("************************************************************")
			print("ERROR: reaction index (reaction.idx) must be greater than 0")
			print("Probably you have to define idx")
			print("************************************************************")
			sys.exit()

		#keeps into account nuclear multeplicity (if option is enabled)
		nuclearMult = "" #default nuclear multeplicity value
		if useNuclearMult:
			uncurledR = [self.reactants[i].name for i in range(len(self.reactants)) if not(self.curlyR[i])]
			if len(uncurledR) == 2:
				if uncurledR[0] == uncurledR[1]: nuclearMult = "0.5d0*"
			if len(uncurledR) == 3:
				if uncurledR[0] == uncurledR[1] and uncurledR[1] == uncurledR[2]: nuclearMult = "0.16666666667d0*"
				if uncurledR[0] == uncurledR[1] and uncurledR[1] != uncurledR[2]: nuclearMult = "0.5d0*"
				if uncurledR[1] == uncurledR[2] and uncurledR[0] != uncurledR[1]: nuclearMult = "0.5d0*"
		self.nuclearMult = nuclearMult
		self.RHS = nuclearMult+"k("+str(self.idx)+")"
		self.RHSvar = "kflux"+str(self.idx)
		ns = []
		actualReactants = self.reactants
		actualCurlyR = self.curlyR
		#if(self.isAutoRev):
		#	actualReactants = self.products
		#	actualCurlyR = self.curlyP
		if useNuclearMult and self.isAutoRev:
			sys.exit("ERROR: nuclear multiplier on reverse reaction is not allowed!")
		i = 0
		for r in actualReactants:
			if len(actualCurlyR) > 0:
				if actualCurlyR[i]: continue #skip curly reactants
			i += 1
			ns.append("n("+str(r.fidx)+")")
			if r.idx <= 0:
				print("************************************************************")
				print("ERROR: species index (molec.idx) must be greater than 0")
				print("Probably you have to define idx")
				print("************************************************************")
				sys.exit()
		if len(ns) > 0: self.RHS += "*" + ("*".join(ns))

	#method:build pseudo hash (for unique reactions to avoid duplicates)
	def build_pseudo_hash(self):
		rname = []
		pname = []
		for r in self.reactants:
			rname.append(r.name)

		for p in self.products:
			pname.append(p.name)
		self.pseudo_hash = ("_".join(sorted(rname)))+"|"+("_".join(sorted(pname)))

	#*********************
	#method: check reaction (mass and charge conservation)
	def check(self,mode="ALL"):
		mass_reactants = mass_products = 0.e0
		charge_reactants = charge_products = 0
		for r in self.reactants:
			mass_reactants += r.mass
			charge_reactants += r.charge

		for p in self.products:
			mass_products += p.mass
			charge_products += p.charge
		if mass_reactants != 0 and (mode == "ALL" or "MASS" in mode):
			if abs(1e0 - mass_products / mass_reactants) > 1e-6:
				print("************************************************")
				print("WARNING: problem with mass conservation in reaction", self.idx)
				print("reaction:",self.verbatim)
				print("reactants:", [r.name for r in self.reactants], mass_reactants)
				print("products:", [p.name for p in self.products], mass_products)
				print("mass ratio (prods/reacts):",mass_products/mass_reactants, "(should be 1.0)")
				print("You can remove this check with the -nomassCheck option")
				print("************************************************")
				a = input("Any key to continue q to quit... ")
				if a == "q": sys.exit()
		if abs(charge_products - charge_reactants)!=0 and (mode=="ALL" or "CHARGE" in mode):
			print("************************************************")
			print("WARNING: problem with charge conservation in reaction", self.idx)
			print("reaction:",self.verbatim)
			print("reactants:", [r.name for r in self.reactants], charge_reactants)
			print("products:", [p.name for p in self.products], charge_products)
			print("You can remove this check with the -nochargeCheck option")
			print("************************************************")
			a = keyb_input("Any key to continue q to quit... ")
			if a == "q":
				sys.exit()

	#calcluate enthalpy of formation
	def enthalpy(self):
		if "krome_kph" in self.krate: return

		if len(self.reactants) == 0 and len(self.products) == 0:
			self.dH = 0e0
			return
		if len(self.reactants) == 0 or len(self.products) == 0:
			print("ERROR: you have called enthalpy calculation")
			print(" with empty reactants and/or products!")
			sys.exit()
		reag = self.reactants  # copy reactants
		prod = self.products  # copy products
		available = True #flag for species availablity in the enthalpy dictionary
		rH = pH = 0.e0 #init reactants and prodcuts enthalpy eV
		#loop on reactants
		for xr in reag:
			xr.name = xr.name.upper()
			#print xr.name
			#if(not(xr.name in deltaH)):
			#	available = False
			#	break
			rH += xr.enthalpy #deltaH[xr.name]

		#loop on products
		for xp in prod:
			xp.name = xp.name.upper()
			#print xp.name
			#if(not(xp.name in deltaH)):
			#	available = False
			#	break
			pH += xp.enthalpy
		self.dH = None
		if available:
			self.dH = (rH-pH)*1.60217657e-12 #eV->erg (cooling<0)


	#*********************
	#alias for enthalpy calculation
	def computeEnthalpy(self):
		self.enthalpy()


	#*********************
	#calculate reverse reaction using polynomials
	def doReverse(self):
		pidx = "(/"+(",".join([x.fidx for x in self.products]))+"/)"
		ridx = "(/"+(",".join([x.fidx for x in self.reactants]))+"/)"
		ndif = len(self.reactants)-len(self.products)
		kk = "("+self.krate+") * revKc(Tgas,"+ridx+","+pidx+")"
		if ndif != 0: kk +=" * (1.3806488d-22 * Tgas)**("+str(ndif)+")"
		return kk

	#***********************
	#get rate F90
	def getRateF90(self, context, varname="k"):
		sTlimit = ""
		hasTlim = (self.hasTlimitMin or self.hasTlimitMax) #Tmin or Tmax are present
		Tlimfound = False #flag to check if endif is needed after the reaction rate
		if self.kphrate is None and context.useTlimits and hasTlim:
			Tlimfound = True #need to close the if statement opened here
			sTlimit = "if("
			if self.hasTlimitMin: sTlimit += "Tgas."+self.TminOp+"."+self.Tmin #Tmin is present
			if self.hasTlimitMin and self.hasTlimitMax: sTlimit += " .and. " #Tmin and Tmax are present
			if self.hasTlimitMax: sTlimit += "Tgas."+self.TmaxOp+"."+self.Tmax #Tmax is present
			sTlimit += ") then\n"
		kstr = "!" + self.verbatim+"\n" #reaction header
		kstr += "\t" + sTlimit + self.ifrate + " "+varname+"("+str(self.idx)+") = " + self.krate #limit+extraif+rate
		if Tlimfound: kstr += "\nend if" #close the if statement for temperature
		return truncF90(kstr, 60,"*") #truncates long reaction rates

#############################################
#find a species with the given name and return it as an object
def searchSpeciesByName(speciesList, speciesName):
	for xsp in speciesList:
		if xsp.name == speciesName: return xsp

	sys.exit("ERROR: species "+speciesName+" not found by searchSpeciesByName!")

###############################
#convert a LEIDEN datafile to KROME format in the build folder
# see http://home.strw.leidenuniv.nl/~ewine/photo/
def LEIDEN2KROME(build_folder,reactant,products):
	try:
		from scipy import interpolate
	except:
		print("ERROR: scipy not installed!")
		print(" This module is necessary to use LEIDEN xsecs.")
		sys.exit()

	prods = sorted([x.name for x in products])
	data_folder = "data/database/leiden_xsecs/"
	fname1 = data_folder+reactant.name+"__"+("_".join(prods))+".dat"
	fname2 = data_folder+reactant.name+"__"+("_".join(prods[::-1]))+".dat"

	if os.path.exists(fname1):
		fname = fname1
	elif os.path.exists(fname2):
		fname = fname2
	else:
		print("")
		print("ERROR: file "+fname1+" OR")
		print(" file "+fname2+" don't exist!")
		print("Species "+reactant.name+" doesn't have a LEIDEN file.")
		print("Search on http://home.strw.leidenuniv.nl/~ewine/photo/index.php?file=pd.php")
		print(" and copy to "+data_folder)
		print("The file should be in the format R__P_P.dat, ")
		print(" where R is the reactant name and P the products,")
		print(" e.g. " + fname1)
		sys.exit()

	print("automatic reaction from LEIDEN database found: " + reactant.name
		  + " -> "+(" + ".join(prods)))

	#columns references: wavelength, absorption, dissociation, and ionization
	header = ["wlen","abs","diss","ion"]
	data = {k: [] for k in header}
	#read data from the file and store the data and the header
	fleiden = open(fname)
	for row in fleiden:
		srow = row.strip()
		if srow == "": continue
		if srow.startswith("#"): continue
		arow = [float(x) for x in srow.split(" ") if(x!="")]
		#append data to array
		for i in range(len(header)):
			data[header[i]].append(arow[i])
	fleiden.close()

	#some constants to convert nm to eV
	clight = 2.99792458e10 #cm/s
	hplanck = 4.135667516e-15 #eV*s

	#nm->cm->eV
	energyList = [clight*hplanck/(x*1e-7) for x in data["wlen"][::-1]]

	#determine if photoionization or not (check if electron in products)
	process = "diss"
	if "E" in prods: process = "ion"
	xsecList = data[process][::-1]

	#create build folder if not exists
	#(this block is also in prepareBuild and dumpNetwork methods
        #from kromeobj.py
	if not os.path.exists(build_folder):
		os.mkdir(build_folder)
		print("Created " + build_folder)
	fout = open(build_folder+"leiden_"+reactant.name+"__"+("_".join(prods))+".org","w")
	fout.write("#original xsecs (cm2) from Leiden database\n")
	fout.write("#http://home.strw.leidenuniv.nl/~ewine/photo/\n")
	for i in range(len(energyList)):
		fout.write(str(energyList[i])+" "+str(xsecList[i])+"\n")
	fout.close()

	#create interpolated function from Leiden file
	finterp = interpolate.interp1d(energyList, xsecList,kind='linear')

	#write the file to the build folder
	foutx = open(build_folder+"leiden_"+reactant.name+"__"+("_".join(prods))+".dat","w")
	imax = 5000 #number of interpolated points
	emax = min(max(energyList), 3e1) #this is the maximum energy limit for creating tables (eV)
	emin = min(energyList)
	foutx.write("#interpolated xsecs, "+str(imax)+" points in range ["+str(emin)+","+str(emax)+"] eV\n")
	#write data to file using a regular xenergy grid
	for i in range(imax):
		xenergy = i*(emax-emin)/(imax-1)+emin
		#if not line-based interpolates with numpy
		xsec = finterp(xenergy)
		foutx.write(str(xenergy)+" "+str(xsec)+"\n")
	foutx.close()


###############################
#convert a SWRI datafile to KROME format in the build folder
def SWRI2KROME(build_folder,reactant,products, Eth):
	import copy
	try:
		from scipy import interpolate
	except:
		print("ERROR: scipy not installed!")
		print(" This module is necessary to use SWRI xsecs.")
		sys.exit()

	prods = [x.name for x in products]
	data_folder = "data/database/swri_xsecs/"
	fname = data_folder+reactant.name+".dat"
	if not os.path.exists(fname):
		print("")
		print("ERROR: file "+fname+" doesn't exist!")
		print("Species "+reactant.name+" doesn't have a SWRI file.")
		print("Search on http://phidrates.space.swri.edu")
		print(" and copy to "+data_folder+"/"+reactant.name+".dat")
		sys.exit()
        # 
	#create build folder if not exists
	#(this block is also in prepareBuild and dumpNetwork methods
        #from kromeobj.py
	if not os.path.exists(build_folder):
		os.mkdir(build_folder)
		print("Created " + build_folder)
	#read data from the file and store the data and the header
	fswri = open(fname)
	okrow = []
	for row in fswri:
		srow = row.strip()
		if srow == "": continue
		arow = [x.strip() for x in srow.split(" ") if x!=""]
		if arow[0] == "Lambda":
			lastLambda = arow[:]
			okrow = []
			continue
		okrow.append(arow)
	fswri.close()

	#replace exited states to ground, e.g. C1D->C
	lold = lastLambda[:]
	for state in ["1D","1S","3P"]:
		lastLambda = [x.replace(state, "") for x in lastLambda]

	#check columns with same products (after state replacing)
	match = [] #index matching
	founds = []
	for i in range(len(lastLambda)-1):
		match.append([])
		for j in range(i+1,len(lastLambda)):
			if lastLambda[i] == lastLambda[j] and j not in founds:
				match[i].append(j)
				founds.append(j)

	#merge (sum) columns with same products (after state replacing)
	removecol = []
	for j in range(len(match)):
		columns = match[j]
		if len(columns) == 0: continue #skip empty matches
		for icol in columns:
			removecol.append(icol)
			for i in range(len(okrow)):
				okrow[i][j] = float(okrow[i][j])+float(okrow[i][icol])

	#write branches slash-separated and include electrons if needed
	lastLambda = [sorted((x.replace("+", "+/E/")).split("/")) for x in lastLambda]
	#remove empty products if any
	lastLambda = [[y for y in x if y!=""] for x in lastLambda]
	lastLambda = [sorted([x.upper() if x=="e" else x for x in block]) for block in lastLambda]

	#search the column that contains the given branch
	psort = sorted(prods)
	if psort not in lastLambda:
		print("ERROR: products "+(" ".join(psort))+" aren't in the SWRI file "
			+ data_folder+reactant.name + ".dat")
		print("Available:", lastLambda[2:])
		print("Looking for:", psort)
		sys.exit()
	data_col = lastLambda.index(psort) #first column matching
	if data_col in removecol: sys.exit("ERROR: problem with SWRI column merging!")

	print("automatic reaction from SWRI database found: "+reactant.name+" -> "+(" + ".join(psort)))

	#some constants to convert AA to eV
	clight = 2.99792458e10  # cm/s
	hplanck = 4.135667516e-15  # eV*s

	#prepare data for interpolation and dump original data for comparison
	xdata = []
	ydata = []
	foutorg = open(build_folder+"swri_"+reactant.name+"__"+("_".join(prods))+".org","w")
	Eth = 1e99 #default threshold energy
	#for loop reads reversed array to have increasing energy
	xsec_old = 0e0
	xmin = xmax = None
	#these 'Dirac' statements are employed to approximate line width (see comments below)
	isDirac = False #boolean to determine if it is reading a line or not
	hasDirac = False #boolean to determine if lines are found or not
	xDirac = [] #list of the lines found
	dDirac = {"freqL": 0e0, "freqR": 0e0, "xsec": 0e0} #template for line (boundaries + xsec)
	#row format is: [E(eV), [xsec(cm2) for each branch]]
	for row in okrow[::-1]:
		if float(row[0]) <= 0e0: continue
		xenergy = clight*hplanck/(float(row[0])*1e-8) #AA->eV
		xsec = float(row[data_col]) #cross section cm2
		#SWRI notation: when xsec is 1e-35 start to read a line
		# with 1e-35 xsec 1e-35 to get the linewidth
		if xsec == 1e-35:
			hasDirac = True
			isDirac = not isDirac
			if xsec_old == 1e-35: isDirac = True #double 1e-35 to open lines block
			#store the left and the right bound of the line
			if isDirac:
				dDirac["freqL"] = xenergy
			else:
				dDirac["freqR"] = xenergy
				xDirac.append(copy.copy(dDirac))
		#store the xsec value when line value (Dirac) is open
		if xsec != 1e-35 and isDirac:
			dDirac["xsec"] = xsec
		xdata.append(xenergy)
		ydata.append(xsec)
		if xsec > 1e-40 and xmin is None: xmin = xenergy
		if xsec > 1e-40: xmax = xenergy
		#find threshold for photoionization
		if xsec == 0e0 and xsec_old != 0e0:
			if "E" in products: Eth = xenergy_old
		foutorg.write(str(xenergy)+" "+str(xsec)+"\n")
		xsec_old = xsec
		xenergy_old = xenergy
	foutorg.close()
	#note: swri uses the sequence "1e-35 value 1e-35" in the xsec
	# to indicate the with of the line and the averaged xsec.
	# the boolean isDirac is true when reading line, hence
	# if not closed (i.e. still False at EOF) rises an error.
	# False because of double 1e-35 to close the line part
	if not isDirac and hasDirac:
		print("ERROR: in SWRI read line with opened token")
		print(" but never closed!")
		print(" check:"+data_folder+reactant.name + ".dat")
		sys.exit()
	#create interpolated function from SWRI file
	fdata = interpolate.interp1d(xdata, ydata,kind='linear')

	if xmax is None or xmin is None or xmax <= xmin:
		print("ERROR: problem with loading cross section for")
		print(reactant.name+" -> "+(" + ".join(psort)))
		print("xmin:", xmin, "xmax:", xmax)
		sys.exit()

	#write the file to the build folder
	foutx = open(build_folder+"swri_"+reactant.name+"__"+("_".join(prods))+".dat","w")
	imax = 5000 #number of interpolated points
	emax = min(xmax, 3e1) #this is the maximum energy limit for creating tables (eV)
	emin = xmin
	#write data to file using a regular xenergy grid
	for i in range(imax):
		xenergy = i*(emax-emin)/(imax-1)+emin
		#if not line-based interpolates with numpy
		if len(xDirac) == 0:
			foutx.write(str(xenergy)+" "+str(fdata(xenergy))+"\n")
		else:
			#if line-based lines are rectangles of average xsec
			myxd = None
			#loop on the lines found to determine the xsec @ xenergy
			for xd in xDirac:
				if xenergy >= xd["freqL"] and xenergy <= xd["freqR"]:
					myxd = xd #store the line at the corresponding xenergy
					break
			#xsec is zero outside the line boundaries
			xsec = 0e0
			if myxd is not None: xsec = myxd["xsec"] #otherwise is the averaged xsec
			foutx.write(str(xenergy)+" "+str(xsec)+"\n")

	foutx.close()

###############################
# NOTE: Using JANAF tables is discouraged because it is unclear on the exact
# calculations/data and it is incompatible with other databases like NASA
# polynomial ones.
# #convert a JANAF datafile to KROME format in the build folder
# def janaf2krome(build_folder, species):
# 	try:
# 		from scipy import interpolate
# 	except:
# 		print "ERROR: scipy not installed!"
# 		print "This module is necessary to use JANAF thermo tables."
# 		sys.exit()
#
# 	name = species.name
# 	filepath = "data/database/janaf/" + name + ".dat"
#
# 	temperature = [] #temperature in K
# 	gibbs = [] #formation gibbs free energy in kJ/mol
# 	# add more lists if needed from janaf file
#
# 	with open(filepath, "r") as janaffile:
# 		for row in janaffile:
# 			srow = row.strip()
# 			if srow == "":
# 				continue
# 			if not is_number(srow[0]):
# 				continue
# 			if "FUGACITY" in srow:
# 				continue
# 			arow = [x.strip() for x in srow.split("\t")]
# 			if len(arow) < 8:
# 				print ("ERROR: there is a problem with row " + str(arow) +
# 				" This is most likely the bug in JANAF tables where spaces are"
# 				"used instead of tabs when going above the fugacity pressure. "
# 				"Please change this in the input table" + filepath
# 				)
# 				sys.exit()
#
# 			temperature.append(float(arow[0]))
# 			# JANAF provides gibbs formation energy w.r.t. reference
# 			# value. To get the uncorreect value, we need to add HminH0.
# 			HminH0 = arow[4] #enthalpy - reference enthalpy in kJ/mol
# 			gib = arow[6] #formation gibbs in kJ/mol
# 			# NOTE: It is unclear from the JANAF database what should be
# 			# used as formation Gibbs free energy (if it needs a correction).
# 			gibbs.append( float(gib) )
# 			# gibbs.append( float(gib) + float(HminH0) )
#
#
# 	#create interpolated function from JANAF file
# 	fdata = interpolate.interp1d(temperature, gibbs, kind='linear')
# 	imax = 200 #number of interpolated points
# 	tmax = max(temperature)
# 	tmin = min(temperature)
# 	#write data to file using a regular xtemp grid
# 	with open(build_folder + "thermo_" + name + ".dat", "w") as fout:
# 		fout.write("#Regularly spaced JANAF table for " + name + "\n")
# 		fout.write("#temperature (K) formation gibbs free energy (kJ/mol)" + "\n")
#
# 		for i in range(imax):
# 			xtemp = i*(tmax-tmin)/(imax-1)+tmin
# 			fout.write(str(xtemp) + " " + str(fdata(xtemp)) + "\n")
#
# 	print "written JANAF " + name + " table"
#

################################
#split molecule name assuming elements written as He
def listA(arg):
	parts = []
	part = ""
	for a in list(arg):
		if a.islower():
			part += a
		else:
			if part != "": parts.append(part)
			part = a
	parts.append(part)
	return parts

####################
#fill with spaces the end of a string up to columns characters
def fillSpaces(instring, columns):
	astring = str(instring)
	if columns < len(astring): return "#"*columns
	return astring+(" "*(columns-len(astring)))

#####################
def compute_Hdata(arg):
	eV2kJmol = 96.4869e0 #eV -> kJ/mol
	#enthalpy data (DH: enthalpy of fomration (kJ/mol), EA: electron affinity (eV)
	#IE: ionization energy (eV) - note: @298K
	#see http://webbook.nist.gov/chemistry/
	xHData = dict()
	xHData["H"] = {"DH":218e0, "IE":13.59844e0, "EA":0.754195e0}
	xHData["H2"] = {"DH":0e0, "IE":15.42593e0}
	xHData["He"] = {"DH":0e0, "IE":24.58741e0}
	xHData["H2+"] = {"DH":1497e0}
	xHData["E"] = {"DH":6.2e0} #5/2RT @298K
	xHData["O"] = {"DH":249.18e0, "IE":13.61806e0, "EA":1.439157e0}
	xHData["O2"] = {"DH":0e0, "IE":12.0697e0, "EA":0.4480e0}
	xHData["OH"] = {"DH":38.99e0, "IE":13.017e0, "EA":1.82767e0}
	xHData["H2O"] = {"DH":-241.826e0, "IE":12.621e0, "EA":12.65e0}
	xHData["H3O+"] = {"DH":603.417e0}
	xHData["H3+"] = {"DH":1.1e3} #from pag.37 Chemistry of the Elements (N. N. Greenwood,A. Earnshaw)
	xHData["C"] = {"DH":716.68e0, "IE":11.26030e0, "EA":1.262114e0}
	xHData["C2"] = {"DH":837.74e0, "IE":11.4e0, "EA":3.273e0}
	xHData["CH"] = {"DH":594.13e0, "IE":10.64e0, "EA":1.26e0}
	xHData["CH2"] = {"DH":386.39e0, "IE":10.396e0, "EA":0.652e0}
	xHData["CH3"] = {"DH":145.59e0, "IE":9.84e0, "EA":0.08e0}
	xHData["HCO"] = {"DH":43.51e0, "IE":8.12e0, "EA":0.313e0}
	xHData["HOC+"] = {"DH":978.7e0} #using reactions 1-2 in Li+2008, J. Chem. Phys.129, 244306
	xHData["CO"] = {"DH":-110.53e0, "IE":14.014e0, "EA":1.32608e0}
	xHData["CO2"] = {"DH":-393.51e0, "IE":13.777e0, "EA":-0.599986e0}
	xHData["N"] = {"DH":472.68e0, "IE":14.534e0}
	xHData["NO"] = {"DH":90.29e0, "IE":9.264e0}
	xHData["CN"] = {"DH":435.14e0, "IE":13.598e0}
	xHData["N2"] = {"DH":0e0, "IE":15.581e0}
	xHData["HCN"] = {"DH":135.15e0, "IE":13.6e0, "EA":0.00156e0}
	xHData["HNC"] = {"DH":207.94} #Wenthold 2000, J. Phys. Chem. A 104, 5612
	xHData["HNO"] = {"DH":99.58e0, "IE":10.1e0, "EA":0.338e0}
	xHData["NH"] = {"DH":376.56e0, "IE":13.49e0, "EA":0.37e0}
	xHData["NH2"] = {"DH":190.37e0, "IE":10.78e0, "EA":0.7710e0}
	xHData["NH3"] = {"DH":-45.9e0, "IE":10.07e0}
	xHData["NH4+"] = {"DH":-132.5e0} #(aq) D.Ebbing, S.D.Gammon, General Chemistry, Enhanced Edition
	xHData["N2H"] = {"DH":251.46e0, "IE":7.8e0} #Cs symmetry, Matus+2006, J. Phys. Chem. A 110, 10116
	xHData["HCNH"] = {"DH":250e0, "IE":9.41e0} #Cowles+1991, J.Chem.Phys. 94, 3517
	xHData["Al"] = {"DH":250e0}
	xHData["AlOH"] = {"DH":-179.91e0}
	xHData["AlO3H3"] = {"DH":-1016.67e0}
	xHData["AlO2H2"] = {"DH":-507.661e0}
	xHData["AlO2H"] = {"DH":-355.472e0}
	xHData["AlO"] = {"DH":42.98e0}
	xHData["AlO2"] = {"DH":-38.658e0}
	xHData["Al2O"] = {"DH":-148.611e0}
	xHData["Al2O2"] = {"DH":-403.096e0}
	xHData["Al2O3"] = {"DH":-546.891e0}


	#extend with uppercase species
	exHData = dict()
	for k,v in xHData.items():
		exHData[k.upper()] = v
	xHData.update(exHData)

	#if DH is present no need to compute it
	if arg in xHData: return xHData[arg]

	#find neutral
	neutral = arg.replace("+","").replace("-","")

	#check if neutral DH exists
	if neutral not in xHData:
		sys.exit("ERROR: cannot compute DH for "+arg+" since "+neutral+" is missing!")

	#do calculation for cation
	if "+" in arg:
		#check if neutral's IE exists
		if "IE" not in xHData[neutral]:
			sys.exit("ERROR: cannot compute DH for "+arg+" since IE "+neutral+" is missing!")
		#DH(cation) = DH(neutral) + ionization energy - electron entalphy
		return {"DH": xHData[neutral]["DH"] + xHData[neutral]["IE"]*eV2kJmol - xHData["E"]["DH"]}

	#do calculation for anion
	if "-" in arg:
		#check if neutral's EA exists
		if "EA" not in xHData[neutral]:
			sys.exit("ERROR: cannot compute DH for "+arg+" since EA "+neutral+" is missing!")
		#DH(anion) = DH(neutral) - electron affinity + electron entalphy
		return {"DH": xHData[neutral]["DH"] - xHData[neutral]["EA"]*eV2kJmol + xHData["E"]["DH"]}

	#if DH doesn't exist for neutral no way to get it
	sys.exit("ERROR: cannot compute DH for "+arg+", it only works for cations/anions")


################################
def compute_DHreact(listRR, listPP):
	DH = 0e0
	return sum([compute_Hdata(x)["DH"] for x in listPP])\
		- sum([compute_Hdata(x)["DH"] for x in listRR])

################################
def generateCustom(readCustomFile):
	from random import random as rand
	from os import listdir
	from os.path import isfile,join
	from math import log10
	import time
	from scipy.optimize import bisect

	#defaults
	custom = dict()
	custom["atoms"] = ["E","H","He"]
	custom["anions"] = "yes"
	custom["cations"] = "yes"
	custom["maxatoms"] = "99"
	custom["maxrea"] = "99"
	custom["maxprod"] = "99"
	custom["photorates"] = "no"
	custom["include"] = []
	custom["exclude"] = []
	custom["present"] = []
	custom["only"] = []

	#list of tokens where an array is expected
	arrayTokens = ["atoms", "include", "exclude", "present", "only"]

	#read custom file and store the info in a dict
	fhcustom = open(readCustomFile)
	for row in fhcustom:
		srow = row.strip()
		if srow.strip() == "": continue
		if srow[0] == "#": continue
		isArray = False
		#search for if this is an array token
		for arr in arrayTokens:
			if arr in srow.split(":"):
				isArray = True
				break
		#read data according to array or value
		if isArray:
			custom[srow.split(":")[0].strip()] = [x.strip() for x in srow.split(":")[1].split(",") if x.strip()!=""]
		else:
			custom[srow.split(":")[0].strip()] = srow.split(":")[1].strip()

	#species available
	amols_org = ["H","H2","H+","H-","He","He+","H2+","H3+","E"]
	amols_org += ["O", "O+", "O-", "O2", "OH", "OH+", "OH-"]
	amols_org += ["H2O", "H2O+", "H3O+"]
	amols_org += ["C", "C+", "C-", "C2", "CH", "CH+", "CH2"]
	amols_org += ["CH2+", "CH3+"]
	amols_org += ["HCO+", "HOC+","HCO","CO","CO+","CO2"]
	amols_org += ["N","NO","CN","N2","HCN","HNC","HNO"]
	amols_org += ["NH","NH2","NH3","N2H+","N2H"]
	amols_org += ["N+","NH+","NH2+","NH3+","NH4+","HCN+","HCNH+"]
	amols_org += ["Al","AlOH","AlO2H2","AlO3H3","AlO2H"]
	amols_org += ["AlO","AlO2","Al2O","Al2O2","Al2O3"]

	#exploded species
	emols_org = ["H","HH","H+","H-","He","He+","HH+","HHH3+","-"]
	emols_org += ["O", "O+", "O-", "O2", "OH", "OH+", "OH-"]
	emols_org += ["HHO", "HHO+", "HHHO+"]
	emols_org += ["C", "C+", "C-", "CC", "CH", "CH+", "CHH"]
	emols_org += ["CHH+", "CHHH+"]
	emols_org += ["HCO+", "HOC+","HCO","CO","CO+","COO"]
	emols_org += ["N","NO","CN","NN","HCN","HNC","HNO"]
	emols_org += ["NH","NHH","NHHH","NNH+","NNH"]
	emols_org += ["N+","NH+","NHH+","NHHH+","NHHHH+","HCN+","HCNH+"]
	emols_org += ["Al","AlOH","AlOOHH","ALOOOHHH","AlOOH"]
	emols_org += ["AlO","AlOO","AlAlO","AlAlOO","AlAlOOO"]

	#convert array present into exploded version
	custom["present"] = [emols_org[amols_org.index(x)] for x in custom["present"]]

	#compute missing DH data using http://webbook.nist.gov/chemistry/ion/#DH prescriptions
	HData = dict()
	for species in amols_org:
		HData[species] = compute_Hdata(species)

	amols = []
	emols = []
	#prepare restricted set according to the options
	for i in range(len(amols_org)):
		#skip anions and cations if reqired
		if ("+" in emols_org[i]) and custom["cations"] == "no": continue
		if ("-" in emols_org[i]) and custom["anions"] == "no": continue
		#skip species if required
		if amols_org[i] in custom["exclude"]: continue
		thisMol = emols_org[i].replace("+","").replace("-","")
		if thisMol == "": thisMol = "E" #fix electron for parsing
		countAtoms = 0 #number of atoms found
		#atoms are sorted for replacing
		for at in sorted(custom["atoms"], key=lambda x:len(x))[::-1]:
			#replace until atom is found in the exploded species
			while at in thisMol:
				countAtoms += 1
				thisMol = thisMol.replace(at,"",1)
		#check max number of atoms per species
		if countAtoms > int(custom["maxatoms"]): continue
		#check residuals from replacing
		if(is_number(thisMol) or thisMol==""):
			amols.append(amols_org[i])
			emols.append(emols_org[i])

	#include additional molecules
	for mol in custom["include"]:
		if mol == "": continue
		if mol in custom["exclude"]: sys.exit("ERROR: "+mol+" present in exclude AND in include on file "+readCustomFile)
		if mol not in amols_org:
			print("ERROR: can't recognize this molecule: "+mol)
			print(" you should add this to the automatic list or correct any typo.")
			sys.exit()
		amols.append(mol)
		emols.append(emols_org[amols_org.index(mol)])

	#check for only species (only reactions with these species)
	if len(custom["only"]) > 0:
		amols = []
		emols = []

	for mol in custom["only"]:
		if mol == "": continue
		if mol not in amols_org:
			print("ERROR: can't recognize this molecule: "+mol)
			print(" you should add this to the automatic list or correct any typo.")
			sys.exit()
		amols.append(mol)
		emols.append(emols_org[amols_org.index(mol)])

	#check number of species found
	if len(amols) == 0:
		print("ERROR: no species available for custom reactions,")
		print(" check options in "+readCustomFile+" file.")
		sys.exit()
	print("Search custom reactions using the following species")
	for i in range(len(amols)//5+1):
		#this simply write 5 species per line
		print(" "+(" ".join(amols[i*5:min((i+1)*5,len(amols))])))

	#temporary file name
	tmpNumber = "network" #str(int(rand()*1e8))
	tmpFname = "cstm"+tmpNumber+".tmp"
	tmpFnameAll = "cstm"+tmpNumber+"_all.tmp"

	print("building automatic product/reactant combinations...")
	#build product/reactant combinations
	combs = []
	for mol1 in emols:
		cc = "".join(sorted(listA(mol1))) #exploded 1body combination
		combs.append([cc,[mol1]]) #single prod/react
		for mol2 in emols:
			if int(custom["maxprod"]) < 2: continue
			if int(custom["maxrea"]) < 2: continue
			cc = "".join(sorted(listA(mol1)+listA(mol2))) #exploded 2body combination
			if "++" in cc: continue #exclude cation-cation reactions
			if "--" in cc: continue #exclude anion-anion reactions
			cc = cc.replace("+-","") #neutralize
			combs.append([cc,sorted([mol1,mol2])]) #double prod/react
			for mol3 in emols:
				if int(custom["maxprod"]) < 3: continue
				if int(custom["maxrea"]) < 3: continue
				#exploded 3body combination
				cc = "".join(sorted(listA(mol1)+listA(mol2)+listA(mol3)))
				#cc = cc.replace("+-","")
				if "+" in cc: continue #exclude cation 3body
				if "-" in cc: continue #exclude anion 3body
				combs.append([cc,sorted([mol1,mol2,mol3])]) #triple prod/react

	#build reactions from product/reactant combinations
	#combinations are [exploded,[species]]
	print("building automatic reactions... (it may take a while)")
	time0 = time.time()
	custRea = []
	lenPresent = len(custom["present"])
	#loop on reactants
	countRR = 0
	for RR in combs:
		if countRR % 5000 == 1:
			print("time to go (s):",round((time.time()-time0)/countRR*(len(combs)-countRR),2))
		countRR += 1
		JRR = ("".join(sorted("".join(RR[1]))))
		lenRR = len(RR[1])
		#loop on products
		for PP in combs:
			#check if exploded are the same
			if RR[0] != PP[0]: continue
			#3body->3body ignored
			if lenRR == 3 and len(PP[1]) == 3: continue
			#anion-cation->cation-anion ignored
			JPP = ("".join(sorted("".join(PP[1]))))
			if "+-" in JRR and "+-" in JPP: continue
			#check if reactant/products are different (avoid A+B->A+B)
			if RR[1] == PP[1]: continue
			#check "present" statement
			if lenPresent > 0:
				allSpecies = RR[1]+PP[1] #all the species R+P
				anyFound = False
				#loop on species
				for sp in allSpecies:
					if sp in custom["present"]:
						anyFound = True
						break
				#if none of the "present" species are found skip reaction
				if not anyFound: continue

			#prepare the reaction
			cRea = sorted([RR,PP])
			#if not already present add
			if cRea not in custRea: custRea.append(cRea)

	print("searching automatic reactions in the dbase...")
	#search reactions in the database
	autoreacts = [] #dbase array contains dictionary with reaction data
	fdbase = "data/database/"
	file_list = [f for f in listdir(fdbase) if isfile(join(fdbase,f))]
	extraVars = dict() #dict of the extra varaibles, with key=filename
	for fname in file_list:
		if "~" in fname: continue #skip temp files
		fhauto = open(fdbase+fname)
		for row in fhauto:
			srow = row.strip()
			if srow.strip() == "": continue
			if srow == "#BREAK DATABASE": break #skip non database files
			if srow[0] == "#": continue
			#skip phtorates if not needed
			if "@photorates:yes" in srow.lower().replace(" ", ""):
				if custom["photorates"] != "yes": break
				continue
			if "@var:" in srow: continue
			if "@type:" in srow: autorea = dict() #begin reaction
			autorea[srow.split(":")[0].replace("@","").strip()] = srow.split(":")[1].strip()
			if "@rate:" in srow: autoreacts.append(autorea) #end reaction

	DHlimit = 3e3 #enthalpy limit to accept a reaction, K
	sDHlimit = str(round(DHlimit/1e1**int(log10(DHlimit)),2))+"d"+str(int(log10(DHlimit)))
	reaFound = []
	fhTmpAll = open(tmpFnameAll,"w")
	fhTmpAll.write("#Reactions automatically found using custom instructions file "+readCustomFile+"\n")
	fhTmpAll.write("#This file also indicates the reactions found in the database (True/False), \n")
	fhTmpAll.write("# satisfying the enthalpy upper limit of "+str(sDHlimit)+" K (*), and that should be checked\n")
	fhTmpAll.write("# because are NOT in the database AND have an enthalpy below that the limit (#).\n")
	fhTmpAll.write(fillSpaces("#IDX",5) + fillSpaces("REACTANS",20) + "    " + fillSpaces("PRODUCTS",20)\
		+ fillSpaces("ENTHALPY (K)",18) + fillSpaces("FOUND IN DBASE",16)\
		+ fillSpaces("<"+str(sDHlimit)+"K",9) + "check\n")

	print("automatic reactions found:",len(custRea)*2)  # including reverse
	print("writing automatic reactions to file... (it may also take a while)")
	time0 = time.time()
	kJmol2K = 120.274e0 #kJ/mol -> K
	#search created reactions in the database
	iCount = cCount = 0
	for crea in custRea:
		if cCount % 500 ==1:
			print("time to go (s):",round((time.time()-time0)/cCount*(len(custRea)-cCount),2))
		cCount += 1
		#prepare custom products/reactants
		CC1 = crea[0][1]
		CC1 = sorted([amols[emols.index(x)] for x in CC1])
		CC2 = crea[1][1]
		CC2 = sorted([amols[emols.index(x)] for x in CC2])
		inDatabaseFwd = inDatabaseRev = False
		DHCC1 = sum([HData[x]["DH"] for x in CC1])
		DHCC2 = sum([HData[x]["DH"] for x in CC2])
		#loop into the databse
		for autorea in autoreacts:
			#prepare database products and reactants
			PP = sorted([x.strip() for x in autorea["prods"].split(",")])
			RR = sorted([x.strip() for x in autorea["reacts"].split(",")])
			#append new reactions found
			if RR==CC1 and PP==CC2 and [CC1,CC2] not in reaFound:
				reaFound.append([CC1,CC2])
				inDatabaseFwd = True
			if RR==CC2 and PP==CC1 and [CC2,CC1] not in reaFound:
				reaFound.append([CC2,CC1])
				inDatabaseRev = True
			if inDatabaseRev and inDatabaseFwd: break
		JJ1 = ("".join(CC1))
		JJ2 = ("".join(CC2))
		#write all the reactions even if not in the database
		RRall = (" + ".join(CC1))
		PPall = (" + ".join(CC2))
		rtype = ""
		if "+" in JJ1 and "-" in JJ1 and len(CC2) == 3: rtype = "+-3prods"
		if list(set(CC1).intersection(CC2)): rtype = "catal"
		if len(CC2) == 1: rtype = "form"
		DH = (DHCC2-DHCC1) * kJmol2K #enthalpy products - reactants
		fav = (" *" if (DH<DHlimit) else "")
		check = (" #" if(DH<DHlimit and not(inDatabaseFwd)) else "")
		if len(CC1) < 3 or rtype == "catal":
			if inDatabaseFwd or DH < DHlimit:
				fhTmpAll.write(fillSpaces(iCount+1,5) + fillSpaces(RRall,20) + " -> " + fillSpaces(PPall,20)\
					+ fillSpaces(DH,18) + fillSpaces(inDatabaseFwd,16)\
					+ fillSpaces(fav,9) + fillSpaces(check,3) + rtype + "\n")
				iCount += 1

		RRall = (" + ".join(CC2))
		PPall = (" + ".join(CC1))
		rtype = ""
		if ("+" in JJ2) and ("-" in JJ2) and (len(CC1) == 3): rtype = "+-3prods"
		if list(set(CC2).intersection(CC1)): rtype = "catal"
		if len(CC1) == 1: rtype = "form"
		DH = (DHCC1-DHCC2) * kJmol2K #enthalpy products - reactants
		fav = (" *" if (DH<DHlimit) else "")
		check = (" #" if(DH<DHlimit and not(inDatabaseRev)) else "")
		if len(CC2) < 3 or rtype == "catal":
			if inDatabaseRev or (DH < DHlimit):
				fhTmpAll.write(fillSpaces(iCount+2,5) + fillSpaces(RRall,20) + " -> " + fillSpaces(PPall,20)\
					+ fillSpaces(DH,18) + fillSpaces(inDatabaseRev,16)\
					+ fillSpaces(fav,9) + fillSpaces(check,3) + rtype + "\n")
				iCount += 1

	fhTmpAll.close()

	#check number of reactions found
	if len(reaFound) == 0:
		print("ERROR: no custom reactions found in the database,")
		print(" check options in "+readCustomFile+" file.")
		sys.exit()

	print("Custom reactions found:",len(reaFound))

	#sort reactions according to their format
	reaFound = sorted(reaFound, key=lambda x:10*len(x[0])+len(x[1]))
	idx = 0
	fmtHashOld = 0
	fhTmp = open(tmpFname, "w")
	print("writing custom reactions on file "+tmpFname)
	fhTmp.write("#this is an automatically-generated network with the options below\n")
	for k,v in custom.items():
		if not isinstance(v, basestring): v = (",".join(v))
		fhTmp.write("#"+k+": "+v+"\n")
	#write results into a network file
	for rea in reaFound:
		fmtHash = 10*len(rea[0])+len(rea[1]) #format hash
		#new format hash requires new format token
		if fmtHashOld != fmtHash: fhTmp.write("\n@format:idx,"+("R,"*len(rea[0]))+("P,"*len(rea[1]))+"rate\n")
		fhTmp.write(str(idx+1)+","+(",".join(rea[0]+rea[1]))+",auto\n")
		fmtHashOld = fmtHash
		idx += 1
	fhTmp.close()

	#this function returns the name of the custom file
	return tmpFname

#################################
#read operators and set attributes
def readTOpt(rea):
	ops = ["GT","LT","GE","LE",">","<"]
	for op in ops:
		if op in rea.Tmin:
			rea.Tmin = rea.Tmin.replace(op,"").replace("..","")
			rea.TminOp = op.replace(">","GT").replace("<","LT")
		if op in rea.Tmax:
			print(rea.Tmax)
			rea.Tmax = rea.Tmax.replace(op,"").replace("..","")
			rea.TmaxOp = op.replace(">","GT").replace("<","LT")
	return rea

#################################
#create tabvar (probably not the best interface ever)
def create_tabvar(mytabvar,mytabpath,mytabxxyy,anytabvars,anytabfiles,anytabpaths,anytabsizes,coevars):
	if mytabvar.split("_")[0].lower() != "user":
		print("ERROR: to avoid conflicts common variables with @tabvar should begin with user_")
		print(" you provided: "+mytabvar)
		print(" it should be: user_"+mytabvar)
		sys.exit()
	#check if file exists
	if not file_exists(mytabpath):
		print("ERROR: file "+mytabpath+" not found!")
		print(" note that the path must be relative to the ./krome command")
		sys.exit()

	#read the size of the table from the first line file
	fhtab = open(mytabpath)
	for tabrow in fhtab:
		stabrow = tabrow.strip()
		if stabrow == "": continue
		if stabrow[0] == "#": continue
		if "," in stabrow:
			mytabsize = [xx.strip() for xx in stabrow.split(",")]
		else:
			print("ERROR: the file "+mytabpath+" must contain the size of the")
			print(" table in the first line (comma separated, e.g. 50,30)")
			sys.exit()
		break
	fhtab.close()

	#retrieve filename from the path
	mytabfile = mytabpath.split("/")[-1] #read the last value

	#store the data in the global arrays
	anytabvars.append(mytabvar)
	anytabfiles.append(mytabfile)
	anytabpaths.append(mytabpath)
	anytabsizes.append(mytabsize)

	anytabx = mytabvar+"_anytabx(:)"
	anytaby = mytabvar+"_anytaby(:)"
	anytabz = mytabvar+"_anytabz(:,:)"
	anytabxmul = mytabvar+"_anytabxmul"
	anytabymul = mytabvar+"_anytabymul"
	tabf = "fit_anytab2D(" + anytabx + ", &\n" + anytaby + ", &\n" + anytabz+", &\n" + anytabxmul\
		   + ", &\n" + anytabymul + ", &\n" + mytabxxyy + ")"
	if mytabvar not in coevars:
		coevars[mytabvar] = [len(coevars), tabf]

	print("Found tabvar:", mytabvar, "("+mytabpath+")", "["+(",".join(mytabsize))+"]")

#############################
def addVarCoe(mytabvar,tabf,coevars):
	if mytabvar not in coevars:
		coevars[mytabvar] = [len(coevars),tabf]

############################
def get_cooling_dict():
	#the keys of this list must be lowercase.
	#the number is the corresponding integer index for the given cooling
	idxcoo = {"H2": 1, "H2GP": 2, "atomic": 3, "CEN": 3, "HD": 4, "Z": 5, "metal": 5, "dH": 6,
			  "enthalpic": 6, "dust": 7, "compton":8, "CIE": 9, "continuum": 10, "cont": 10,
			  "exp": 11, "expansion": 11, "ff": 12, "bss": 12, "custom": 13, "CO":14,
			  "ZCIE": 15, "ZCIENOUV": 16, "ZExtend": 17, "GH": 18, "OH": 19, "H2O": 20, "HCN": 21}
	idxcoo = {k.lower(): v for k, v in idxcoo.items()}
	return idxcoo

#############################
#cooling index list
def get_cooling_index_list():

	idxcoo = get_cooling_dict()
	#loop on the index to write variables as idx_cool_H2 = 1
	idxscoo = []
	maxv = 0 #maximum index found is the size of the cooling array
	for k, v in idxcoo.items():
		idxscoo.append([v, "idx_cool_"+k+" = "+str(v)])
		maxv = max(maxv, v)
	idxscoo = sorted(idxscoo, key=lambda x: x[0])
	idxscoo.append([99, "ncools = "+str(maxv)])
	return [x[1] for x in idxscoo]

############################
def get_heating_dict():
	#the number is the corresponding integer index for the given heating
	idxhea = {"chem": 1, "compress": 2, "compr": 2, "photo": 3, "dH": 4, "enthalpic": 4,
			  "photoAv": 5, "Av": 5, "CR": 6, "dust": 7, "xray": 8, "visc": 9,"viscous": 9,
			  "custom": 10, "ZCIE": 11}
	idxhea = {k.lower(): v for k, v in idxhea.items()}
	return idxhea

#############################
#heating index list
def get_heating_index_list():

	idxhea = get_heating_dict()
	idxshea = []
	maxv = 0
	for k, v in idxhea.items():
		idxshea.append([v, "idx_heat_"+k+" = "+str(v)])
		maxv = max(maxv, v)
	idxshea = sorted(idxshea, key=lambda x: x[0])
	idxshea.append([99, "nheats = "+str(maxv)])
	return [x[1] for x in idxshea]


####################################
#solar metallicities
def get_solar_abundances(fileName="data/asplund.dat", keyName="Solar"):
	#solar abundances from Tab.1 in Asplund+2009
	# following their definition
	#  log10(epsilon) = log10(n/nH) + 12
	# where n is the number densiity of the given element,
	# while nH is the number density of H
	# This function returns log10(n/nH)

	# Z:     Atomic number [I3]
	# A:     Average weight [F8.3]
	# Na:    Element name [A3]
	# Solar: Epsilon abundance in the solar photosphere [F6.2]
	# Err:   Error on epsilon value [F6.2]
	# Meteo: Epsilon abundance in meteorites (specifically CI chondrites) [F6.2]
	# Err:   Error on epsilon value [F6.2]

	#format spacing and names
	fmtSpace = [3, 8, 3] + [6]*4
	fmtName = ["Z", "A", "Name", "Solar", "ErrSolar", "Meteor", "MeteorErr"]

	solar_abs = dict()

	fhz = open(fileName)
	#loop on data file
	for row in fhz:
		srow = row.strip()
		#skip comments
		if srow == "": continue
		if srow.startswith("#"): continue
		#fill trailing with large number of spaces to avoid format problems
		srow = srow+(""*sum(fmtSpace))
		#cursor position
		istep = 0
		data = dict()
		#loop on format pieces
		for i in range(len(fmtSpace)):
			#get size of data
			nsp = fmtSpace[i]
			#store data in dict
			data[fmtName[i]] = srow[istep:istep+nsp]
			#increase cursor position
			istep += nsp
		#if data are not empty store abundance
		if data[keyName].strip() != "":
			solar_abs[data["Name"].strip()] = float(data[keyName])

	if keyName == "Solar":
		solar_out = dict()
		for k, v in solar_abs.items():
			solar_out[k] = str(v-12e0)
		return solar_out
	else:
		return solar_abs

###############################
#get abundances assuming specific depletion and dust/gas mass ratio
def getAbundancesWithDepletion(d2g = .01):

	#species not depleted (fdep=0)
	nodep = ["H","He","Ne","Ar","Kr","Xe","Rn","Hg"]
	#species partially depleted (0<fdep<1)
	dep = {"O":0.8,"N":0.5,"Si":0.8}
	#species to compute depletion
	depFree = "C"

	#list of species with fdep=0
	nodepList = nodep+[depFree]

	#get Asplund exponents of solar number densities, i.e. log(epsilon)-12
	Zaspl = get_solar_abundances()
	#get average atomic mass
	Waspl = get_solar_abundances(keyName="A")
	#compute Z in mass for each atom (not normalized)
	Zmassn = {k:(Waspl[k]*1e1**float(v)) for (k,v) in Zaspl.items()}
	#total Z not normalized
	Ztot = sum(Zmassn.values())
	#normalize to have Z
	Zmass = {k:v/Ztot for (k,v) in Zmassn.items()}

	#create depletion factors dictionary (1=all in dust)
	fdep = {k:0e0 for k in nodepList}
	fdep.update({k:1e0 for k in Zmass.keys() if(not(k in nodepList))})
	fdep.update(dep)

	#compute X,Y,Z
	X = Zmass["H"]
	Y = Zmass["He"]
	Z = sum([v for (k,v) in Zmass.items() if(not(k in ["H","He"]))])

	#Z depleted (excluded unknown)
	Zdep = sum([v*fdep[k] for (k,v) in Zmass.items()])

	#Z in dust
	Zd = d2g/(d2g+1e0)
	#Z in gas
	Zg = Z - Zd

	#Z of partially depleted species
	Zx = Zd-Zdep

	#get depletion of unknown species
	ff = Zx/Zmass[depFree]

	#trigger error if required depletion > 1 or negative
	if ff > 1e0 or ff < 0e0:
		print("X:",X,"Y:",Y, "Z:", Z)
		print("Zd:",Zd, "Zg:",Zg, "Zd/Z:",Zd/Z, "Zg/Z:",Zg/Z)
		print(dep)
		print(depFree+":", ff)
		sys.exit("ERROR: problems with "+depFree+" depletion!")
	fdep[depFree] = ff

	nH = Zmass["H"]*(1.-fdep["H"])/Waspl["H"]
	return {k:v*(1.-fdep[k])/Waspl[k]/nH for (k,v) in Zmass.items()}
	#print {k:v/Waspl[k] for (k,v) in Zmassn.iteritems()}


#*****************************
def get_Ebind(fileName="data/Ebare_ice.dat",surface="bare"):
	fh = open(fileName)

	Ebind = dict()
	for row in fh:
		srow = row.strip()
		if srow.startswith("#"): continue
		if srow == "": continue
		(name, Ebare, Eice) = [x for x in srow.split(" ") if(x!="")]
		if surface == "bare":
			Ebind[name] = float(Ebare)
		elif surface == "ice":
			Ebind[name] = float(Eice)
		else:
			sys.exit("ERROR: unknown surface type "+surface)

	return Ebind

###################################
#vibrational constant dictionary
#IRIKURA J. Phys. Chem. Ref. Data, Vol. 36, No. 2, 2007
#energy in cm-1, returns K
#returns False if arg not found in list
def get_ve_vib(arg):
	ve = {"H2":4401.213,
		"HD":3813.15,
		"D2":3115.5,
		"C2":1855.0663,
		"C2-":1781.189,
		"CH":2860.7508,
		"CO":2169.75589,
		"CO+":2214.127,
		"N2":2358.57,
		"N2+":2207.0115,
		"NH":3282.72,
		"NO":1904.1346,
		"NO+":2376.72,
		"O2":1580.161,
		"O2+":1905.892,
		"OH":3737.761}
	if arg in ve:
		return ve[arg]*1.42879e0 #cm-1 to K
	else:
		return None

###################################
#rotational constant Be dictionary
#from NIST and Atkins Book
#constant in cm-1, returns K
#returns False if arg not found in list
def get_be_rot(arg):
	ve = {"H2":60.853,
		"H2+":42.9,
		"HD":45.644,
		"D2":30.443,
		"N2":1.9982,
		"O2":1.4264,
		"CO":2.78}
	if arg in ve:
		return ve[arg]*1.42879e0 #cm-1 to K
	else:
		return None

# Equivalent radius of a species in cm
def get_radius(arg):
	radius = {"TIO2": 1.62e-8,
	# Interatomic distance from Jeong et al 2000 DOI:10.1088/0953-4075/33/17/319
			"AL2O3": 3.304e-8,
	# Interatomic distance O-Al-O (linear geometry) from Archibong et al 1999
	# doi: 10.1021/jp983695n
			"MGO":  0.865e-8,
	# Half a bond length form Farrow et al 2014 doi:10.1039/C4CP01825G
			"SIO": 0.75765e-8
	# Bond length from Bromley et al 2016 doi:10.1039/c6cp03629e
			}
	if arg in radius:
		return radius[arg]
	else:
		return None

##################################
#check if a file exists
def file_exists(fname):
	import os
	return os.path.exists(fname)

####################################
def int_to_roman(input):
	#Convert an integer to Roman numerals.
	#from http://code.activestate.com/recipes/81611-roman-numerals/
	if type(input) != type(1):
		print("ERROR: expected integer, got " + type(input))
		sys.exit()
	if not 0 < input < 4000:
		print("ERROR: Argument must be between 1 and 3999")
		sys.exit()
	ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
	nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
	result = ""
	for i in range(len(ints)):
		count = int(input / ints[i])
		result += nums[i] * count
		input -= ints[i] * count
	return result

##################################
#store file in a list
def store_file(fname):
	fh = open(fname,"rb")
	fle = []
	for x in fh:
		fle.append(x)
	fh.close()
	return fle

##################################
#restore file from a list
def restore_file(fname,fle):
	fout = open(fname,"w")
	for x in fle:
		fout.write(x)
	fout.close()

##################################
def get_terminal_size(fd=1):
	"""
    Returns height and width of current terminal. First tries to get
    size via termios.TIOCGWINSZ, then from environment. Defaults to 25
    lines x 80 columns if both methods fail.

    :param fd: file descriptor (default: 1=stdout)
	from bit.ly/HteEcQ
    """
	try:
		import fcntl, termios, struct
		hw = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
	except:
		try:
			hw = (os.environ['LINES'], os.environ['COLUMNS'])
		except:
			hw = (25, 80)

	return hw

##################################
#return an example for test.f90
def get_example(nsp, useX):
	sfile = """
		!###################################################
		! WARNING:This is a test auto-generated by KROME, in order to
		! show a bare-minimal code to call the KROME's subroutine.
		! Most of the values could not be appropriate for your
		! problem, since this test is only intended as a general
		! purpose example.

		program test
			use krome_main !use krome (mandatory)
			use krome_user !use utility (for krome_idx_* constants and others)
			implicit none
			integer,parameter::nsp=krome_nmols !number of species (common)
			real*8::Tgas,dt,x(nsp)@rho@,spy

			spy = 3.65d2 * 2.4d1 * 3.6d3 !seconds per year

			call krome_init() !init krome (mandatory)

			x(:) = 1d-20 !default abundances
			x(krome_idx_H) = @nH@ !hydrogen initial abundance
			@norm@

			Tgas = 1d3 !gas temperature (K)
			dt = 1d6 * spy !time-step (s)
			@rho_init@

			!call the solver
			call krome(x(:),@rhof@ Tgas, dt) !call KROME

			print *,"Test OK!"

		end program test

	"""
	if useX:
		sfile = sfile.replace("@rho@",",rho").replace("@nH@","1.d0").replace("@norm@","x(:) = x(:) / sum(x) !normalize")
		sfile = sfile.replace("@rho_init@","rho = 1d-18 !gas density (g/cm3)").replace("@rhof@"," rho,")
	else:
		sfile = sfile.replace("@rho@","").replace("@nH@","1.d4").replace("@norm@","")
		sfile = sfile.replace("@rho_init@","").replace("@rhof@","")
	sfile = sfile.replace("@nsp@",str(nsp))

	if "@" in sfile:
		print(sfile)
		print("ERROR: missing replacement in get_example() function!")
		sys.exit()
	return sfile


#########################################
#parse for variables in a f90 expression
def parsevar(arg):
	arg = arg.lower()
	tks = ["+","-","*","/","(",")"]
	for tk in tks:
		arg = arg.replace(tk,"@")
	while "@@" in arg:
		arg = arg.replace("@@","@")
	#tks2 = ["exp","sqrt","log","log10"]
	return arg.split("@")

##################################
#extend the list slist with the temperature shortcuts
# for the reaction rea
def get_Tshortcut(rea, slist, cvars=[]):
	shcut = ["logT = log10(Tgas) !log10 of Tgas (#)",
	"lnT = log(Tgas) !ln of Tgas (#)",
	"Te = Tgas*8.617343d-5 !Tgas in eV (eV)",
	"lnTe = log(Te) !ln of Te (#)",
	"T32 = Tgas*0.0033333333333333335 !Tgas/(300 K) (#)",
	"t3 = T32 !alias for T32 (#)",
	"t4 = Tgas*1d-4 !Tgas/1d4 (#)",
	"invT = 1.d0/Tgas !inverse of T (1/K)",
	"invT32 = 1.d0/T32 !inverse of T32 (1/K)",
	"invTgas = 1.d0/Tgas !inverse of T (1/K)",
	"invTe = 1.d0/Te !inverse of T (1/eV)",
	"sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)",
	"invsqrT32 = 1.d0/sqrt(T32)",
	"invsqrT = 1.d0/sqrTgas",
	"sqrT32 = sqrt(T32)",
	"Tgas2 = Tgas*Tgas",
	"Tgas3 = Tgas2*Tgas",
	"Tgas4 = Tgas3*Tgas",
	"T0 = 288d0 !standard temperature (K)",
	"T02 = T0*T0",
	"T03 = T02*T0",
	"T04 = T03*T0",
	"T0inv = 1.d0/T0"]

	#split the first part of the shortcut and uses it as a key (e.g. t = Tgas -> t)
	sckey = [(x.split("="))[0].strip().lower() for x in shcut]
	shcut = sorted(shcut, key=lambda x:len(x.split("=")[0].strip()), reverse=True)

	#loop on the shortcuts to find if the rate coefficient employs them
	krea = rea.krate.replace(" ","")
	for x in shcut:
		ax = x.split("=") #split the shortcut
		xvar = parsevar(krea) #parse the variable in the rate coefficient
		krea = krea.replace(ax[0].strip(),"")
		if ax[0].strip().lower() in xvar and x not in slist:
			#search for dependencies between shortcuts
			xtmp = x
			for xx in shcut:
				axx = xx.split("=") #split the shortcut
				xxvar = parsevar(xtmp)  #parse the variables in the rate coefficient found
				xtmp = xtmp.replace(axx[0].strip(),"")
				if axx[0].strip().lower() in xxvar and xx not in slist:
					slist.append(xx) #append the dependent shortcut
			slist.append(x) #append the main shortcut

	#sort the shortcuts by using the index in the list sckey to follow the hierarchy
	slist = sorted(slist, key=lambda x:sckey.index((x.split("="))[0].strip().lower()))

	#remove shortcuts that are already use as user-defined variables
	slistu = []
	for x in slist:
		xkey = x.split("=")[0] #get the variable name
		xFound = False
		#loop on user-defined variables
		for cv in cvars:
			#when variable is found skip it
			if xkey.lower().strip() == cv.lower().strip():
				xFound = True
				break
		if not xFound:
			slistu.append(x)

	#keep only unique shortcuts (remove duplicates)
	ulist = []
	for x in slistu:
		if x not in ulist: ulist.append(x)

	return ulist

##################################
#get list of available commands
def get_usage():
	print("use -h to see the help")
	sys.exit()

##################################
#truncate F90 expression using sep as separator for blocks shorter than sublen
def truncF90(mystr, sublen, sep):
	#split (&\n) the string mystr in parts smaller than sublen using sep as separator
	if mystr.strip() == "": return mystr
	mystr = mystr.replace("**","##")
	mystr = mystr.replace("(/","###")
	mystr = mystr.replace("/)","####")
	astr = mystr.split(sep)
	s = z = ""
	first = True
	for x in astr:
		if len(z+x) > sublen and not first:
			s += "&\n"
			z = ""
		zep = sep
		if first: zep = ""
		s += zep + x
		z += zep + x
		first = False
	return s.replace("####","/)").replace("###","(/").replace("##","**")


##################################
#look for array definition in var token and prepare variable name according to array size
def coeVarArray(varin):
	if "[" in varin:
		varin = varin.replace(" ","") #replace spaces
		var_array_size = varin.split("[")[1].split("]")[0] #grep inside brackets
		varin = varin.replace("["+var_array_size+"]","") #replace brackets and content
		varin = varin+"("+var_array_size+")" #set variable definition

	return varin


################################
def at_extract(arg):
	aarg = arg.split(":")
	if len(aarg) < 2:
		print("ERROR: @label:value format not respected for")
		print(arg)
		sys.exit()
	aarg[1] = (":".join(aarg[1:]))
	aarg[0] = aarg[0].replace("@","")
	return {aarg[0]:aarg[1]}

##################################
#double format for f90 expressions
def format_double(snum):
	snum = str(snum)
	#format string number to F90 double
	if "d" in snum: return snum
	if "e" in snum: return snum.replace("e","d")
	return snum+"d0"

##################################
#format the subelement according to its size (Fe, N, ...)
def format_subel(subel):
	if len(subel) == 1:
		fsubel = subel.upper()
		if fsubel == "G": fsubel = "g"
	elif(len(subel)==2):
		fsubel = subel[0].upper() + subel[1].lower()
	else:
		fsubel = subel
	return fsubel

##################################
#check if the value s is a number or not (return logical)
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

##################################
#parse molecule name using dictionary and atoms list
def parser(name, mass_dic, atoms, thermo_data, dustIdx=0):

	mymol = molec() #oggetto molec
	namecp = name.upper()
	if namecp == "E-": namecp = "E" #avoid double negative charge
	ename = [] #exploded name
	mass = 0. #init mass
	is_atom = True #atom flag
	founds = 0 #atoms found
	mymol.parentDustBin = dustIdx #this species belongs to this dust bin (1-based)

	#if you change these check the same values in kromeobj
	#(employed here for computing number of neutrons)
	me = 9.10938188e-28 #electron mass (g)
	mp = 1.67262158e-24 #proton mass (g)
	mn = 1.6725e-24 #neutron mass (g)


	#atomic number dictionary (add here atoms if needed)
	zdic = {"E":0,
		"D":1,
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
		"TI":22,
		"CR":24,
		"MN":25,
		"FE":26,
		"NI":28}

	#add isotopes Z numbers (same as non-isotopes)
	zdic_copy = zdic.copy()
	for k,v in zdic_copy.items():
		for i in range(2,56):
			zdic[str(i)+k] = v

	#look for species on grain surface
	if "_dust" in name.lower():
		mymol.is_surface = True

	#check for fake species
	if "FK" in name:
		mymol.name = name #name
		mymol.mass = 0e0 #mass (g)
		mymol.ename = name #exploded name
		mymol.charge = 0 #charge
		mymol.zatom = 0 #atomic number
		mymol.fname = name #f90 name
		mymol.nameLatex = name
		mymol.is_atom = True #atom flag
		mymol.fidx = "idx_"+name #f90 index
		mymol.neutrons = 0 #number of neutrons
		return mymol

	#chemisorbed species
	if "_c_dust" in name.lower():
		mymol.is_chemisorbed = True

	clusterables = ["TIO2", "MGO", "SIO", "AL2O3"]
	if name in clusterables:
		mymol.is_clustarable = True
		mymol.radius = get_radius(name)

	#when belongs to dust+idx remove _dust in the name
	if dustIdx>0 and not mymol.is_surface:
		print("ERROR: in parser, dustIdx>0 with a non-surface species")
		print(dustIdx, name)
		sys.exit()

	zatom = 0 #atomic number init
	#loop over charcters
	for atm in atoms:
		a = atm.upper() #capitalize name
		if a not in namecp: continue #skip
		#loop to up to _30 subscript
		for j in range(30):
			if a in namecp:
				idx = namecp.find(a) #find position
				subs = a
				mult = "0" #multiplicator
				for i in range(idx+len(a),len(namecp)):
					if not is_number(namecp[i]): break
					mult += namecp[i] #find multiplicator
					subs += namecp[i]
				imult = max(int(mult),1) #evaluate multiplicator (must be >0)
				mass += mass_dic[a]*imult #compute mass
				if a in zdic: zatom += zdic[a]*imult #increase atomic number
				if format_subel(a) in mymol.atomcount:
					mymol.atomcount[format_subel(a)] += imult #increase atom count
				else:
					mymol.atomcount[format_subel(a)] = imult #init atom count
				ename += [format_subel(a)]*imult #exploded name
				namecp = namecp.replace(subs,"",1) #remove found in name
				if namecp == "": break #if nothing more to find break loop
		if a != "+" and a != "-": founds += imult #count found atoms for is_atom
	if founds > 1: is_atom = False #atoms have only one atom (viz.)

	mymol.atomcount2 = dict()
	natoms = 0
	for  k, v in mymol.atomcount.items():
		if v > 0:
			mymol.atomcount2[k] = v
			if k != "+" and k != "-": natoms += v

	#print name,mymol.atomcount2

	mymol.natoms = natoms #number of atoms (e.g. diatom=2)

	#get vibrational constant in K
	if get_ve_vib(name):
		mymol.ve_vib = get_ve_vib(name)
	#get rotational constant in K
	if get_be_rot(name):
		mymol.be_rot = get_be_rot(name)

	#when dust index changes the name
	if dustIdx > 0: name = name+"_"+str(dustIdx)

	mymol.name = name #name
	mymol.mass = mass #mass (g)
	mymol.ename = sorted(ename) #exploded name
	mymol.charge = 0 #charge
	mymol.zatom = zatom #atomic number
	mymol.fname = name.replace("+","j").replace("-","k") #f90 name
	#name latex version, replace number with subscript, signs with superscripts
	expName = list(name)
	repName = []
	#loop on name parts
	for part in expName:
		if is_number(part): part = "$_"+part+"$"
		if part == "+" or part == "-": part = "$^"+part+"$"
		repName.append(part)
	mymol.nameLatex = ("".join(repName))
	#replace latex name if electron
	if name == "E": mymol.nameLatex = "e$^-$"
	#cooling name is only for atoms, e.g. CIV
	mymol.coolname = name
	if is_atom:
		mymol.coolname = getRomanName(name)

	mymol.is_atom = is_atom #atom flag
	f90idx = "idx_"+name.replace("+","j").replace("-","k").replace("(","_").replace(")","").replace("[","").replace("]","_") #f90 index
	if f90idx.endswith("_"): f90idx = f90idx[:-1] #remove last underscore if any
	mymol.fidx = f90idx #index in f90 format

	if "+" in name: mymol.charge = name.count("+") #get + charge
	if "-" in name: mymol.charge = -name.count("-") #get - charge

	#number of neutrons (computed using total mass)
	Nn = round((mymol.mass - (me*(mymol.zatom - mymol.charge) + mp*(mymol.zatom))) / mn,0)
	mymol.neutrons = int(Nn)

	#name for photoionization reactions (e.g. Sijjj = Sij3)
	jj = kk = ""
	if mymol.charge == 1: jj = "j"
	if mymol.charge > 1: jj = "j"+str(mymol.charge)
	if mymol.charge == -1: kk = "k"
	if mymol.charge < -1: kk = "k"+str(mymol.charge)
	mymol.phname = name.replace("+","").replace("-","") + jj + kk

	#electron has negative charge
	if mymol.name == "E": mymol.charge = -1

	#thermal data
	if mymol.name in thermo_data:
		if "NASA" in thermo_data[mymol.name]:
			mymol.poly1_nasa = thermo_data[mymol.name]["NASA"][10:] #NASA polynomials lower T interval (min-med)
			mymol.poly2_nasa = thermo_data[mymol.name]["NASA"][3:10] #NASA polynomials upper T interval (med-max)
			mymol.Tpoly_nasa = thermo_data[mymol.name]["NASA"][0:3] #(K) [min,med,max] T interval limits

		if "NIST" in thermo_data[mymol.name]:
			mymol.poly1_nist = thermo_data[mymol.name]["NIST"][2:9] #NIST polynomials lower T interval (min-med)
			mymol.Tpoly_nist = thermo_data[mymol.name]["NIST"][0:2] + [0]  #(K) [min,max,0] T interval limits
			#check for multiple temperature ranges
			if len(thermo_data[mymol.name]["NIST"]) > 9:
				mymol.Tpoly_nist[-1] = thermo_data[mymol.name]["NIST"][10]  #(K) [min,med,max] T interval limits
				mymol.poly2_nist = thermo_data[mymol.name]["NIST"][11:] #NIST polynomials upper T interval

	# Check if species has Gibbs free energy data
	# priority for tables manually added to KROME
	# over some database
	GFE_path = "data/thermochemistry/" + mymol.name + ".gfe"
	# thermoJanafTabPath = "data/database/janaf/" + mymol.name + ".dat"
	if file_exists(GFE_path):
		mymol.has_GFE_table = True
	# NOTE: Using JANAF tables is discouraged because it is unclear on the exact
	# calculations/data and it is incompatible with other databases like NASA
	# polynomial ones.
	# elif file_exists(thermoJanafTabPath):
	# 	mymol.has_GFE_table = True
	# 	mymol.hasJanafThermoTable = True

	#compute enthaly @300K using NASA poly
	if mymol.Tpoly_nasa[1] < 3e2:
		p = mymol.poly1_nasa #copy polynomials in the lower range
	else:
		p = mymol.poly2_nasa #copy poly in the upper range
	Tgas = 300. #K
	polyH = p[0] + p[1]*0.5*Tgas + p[2]*Tgas**2/3. + p[3]*Tgas**3*0.25 + p[4]*Tgas**4*0.2 + p[5]/Tgas
	mymol.enthalpy = polyH*8.314472e-3*Tgas*0.01036410e0 #eV

	#checks parsing results
	if len(namecp) > 0:
		print("************************************************")
		print("ERROR: Parsing problem for", name)
		print("Unknown subelements in substring \"" + namecp +"\".")
		print("Probably you have to add some subelements to the dictionary mass_dic.")
		if len(atoms) < 30:
			print("Dictionary now contains the following subelements:")
			print(atoms)
		print("************************************************")
		sys.exit()
	return mymol

###################################
#convert the name of a species to roman name, e.g. C++ to CIII
# also for anion: C- to CmI
def getRomanName(argmetal):
	#cation
	if "+" in argmetal:
		mname = argmetal.replace("+","") + int_to_roman(argmetal.count("+")+1)
	#anion
	elif "-" in argmetal:
		mname = argmetal.replace("-","") + "m"+int_to_roman(argmetal.count("-")+1)
	#neutral
	else:
		mname = argmetal+"I"
	return mname


######################################
#return a list of the core files of krome to check the consistency
def get_file_list():
	files = []
	files.append("krome")
	files.append("kromeobj.py")
	files.append("kromelib.py")
	files.append("patches")
	files.append("tests")
	files.append("tests/test.f90")
	files.append("tests/Makefile")
	files.append("tests/MakefileF90")
	files.append("build")
	files.append("tools")
	files.append("src")
	files.append("src/krome_constants.f90")
	files.append("src/krome_tabs.f90")
	files.append("src/kromeF90.f90")
	files.append("src/krome_user.f90")
	files.append("src/krome_stars.f90")
	files.append("src/krome_user_commons.f90")
	files.append("src/krome_getphys.f90")
	files.append("src/krome_grfuncs.f90")
	files.append("src/krome_phfuncs.f90")
	files.append("src/krome_gadiab.f90")
	files.append("src/krome_subs.f90")
	files.append("src/krome_dust.f90")
	files.append("src/krome_cooling.f90")
	files.append("src/krome_photo.f90")
	files.append("src/krome_ode.f90")
	files.append("src/krome_reduction.f90")
	files.append("src/krome.f90")
	files.append("src/krome_heating.f90")
	files.append("src/krome_commons.f90")
	files.append("data")
	files.append("solver")
	files.append("solver/nleq_all.f")
	files.append("solver/opkda2.f")
	files.append("solver/opkdmain.f")
	files.append("solver/dvode_f90_license.txt")
	files.append("solver/opkda1.f")
	files.append("networks")
	return files

###################################
# modified from:
# http://akiscode.com/articles/sha-1directoryhash.shtml
# Copyright (c) 2009 Stephen Akiki
# MIT License (Means you can do whatever you want with this)
#  See http://www.opensource.org/licenses/mit-license.php

def GetHashofDirs():
	import hashlib, os
	SHAhash = hashlib.sha1()
	fles = get_file_list()
	for fle in fles:
		try:
			f1 = open(fle, 'rb')
		except:
			continue
		buf = f1.read(4096)
		SHAhash.update(hashlib.sha1(buf).hexdigest())

	return SHAhash.hexdigest()

##################################
def clear_dir(folder):
	import os
	for the_file in os.listdir(folder):
		file_path = os.path.join(folder, the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception as e:
			print(e)


#################################
# returns the license of KROME
def get_licence_header(version, codename, short=False):
	import datetime
	header =  """!!*************************************************************
	!! This file has been generated with:
	!! KROME #version# "#codename#" on #date#
	!! Changeset #changeset#
	!!
	!! KROME is a nice and friendly chemistry package for a wide range of
	!! astrophysical simulations. Given a chemical network (in CSV format)
	!! it automatically generates all the routines needed to solve the kinetic
	!! of the system, modelled as system of coupled Ordinary Differential
	!! Equations.
	!! It provides different options which make it unique and very flexible.
	!! Any suggestions and comments are welcomed. KROME is an open-source
	!! package, GNU-licensed, and any improvements provided by
	!! the users is well accepted. See disclaimer below and GNU License
	!! in gpl-3.0.txt.
	!!
	!! more details in http://kromepackage.org/
	!! also see https://bitbucket.org/krome/krome_stable
	!!
	!! Written and developed by Tommaso Grassi
	!!  tgrassi@nbi.dk,
	!!  Starplan Center, Copenhagen.
	!!  Niels Bohr Institute, Copenhagen.
	!!
	!! and Stefano Bovino
	!!  stefano.bovino@uni-hamburg.de
	!!  Hamburger Sternwarte, Hamburg.
	!!
	!! Contributors:
	!! #contributors#
	!!
	!!
	!! KROME is provided \"as it is\", without any warranty.
	!! The Authors assume no liability for any damages of any kind
	!! (direct or indirect damages, contractual or non-contractual
	!! damages, pecuniary or non-pecuniary damages), directly or
	!! indirectly derived or arising from the correct or incorrect
	!! usage of KROME, in any possible environment, or arising from
	!! the impossibility to use, fully or partially, the software,
	!! or any bug or malefunction.
	!! Such exclusion of liability expressly includes any damages
	!! including the loss of data of any kind (including personal data)
	!!*************************************************************\n"""

	if short: header = """!!*************************************************************
	!! This file has been generated with:
	!! KROME #version# on #date#
	!! Changeset #changeset#
	!! see http://kromepackage.org
	!!
	!! Written and developed by Tommaso Grassi and Stefano Bovino
	!!
	!! Contributors:
	!! #contributors#
	!! KROME is provided \"as it is\", without any warranty.
	!!*************************************************************\n"""

	#list of contributors
	contribs = ["D.Galli", "F.A.Gianturco", "T.Haugboelle", "A.Lupi", "J.Prieto", "J.Ramsey", \
		"D.R.G.Schleicher", "D.Seifried", "E.Simoncini", "E.Tognelli", "T.Frostholm", "J.Boulangier"]

	#sort alphabetically
	contribs = sorted(contribs, key=lambda x:x.split(".")[-1])

	#divide in rows with limited length
	contributors = rowtmp = ""
	for contr in contribs:
		rowtmp += contr+", "
		contributors += contr+", "
		if len(rowtmp) > 60:
			rowtmp = ""
			contributors += "\n\t!! "
	#remove the last space+comma
	contributors = contributors[:-2]

	#name of the git master file
	masterfile = ".git/refs/heads/master"
	changeset = ("x"*7) #default unknown changeset
	#if git master file exists grep the changeset
	if(file_exists(masterfile)):
		changeset = open(masterfile).read()

	# replace comments pragmas
	datenow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	header = header.replace("#date#", datenow).replace("#version#", version)
	header = header.replace("#codename#", codename).replace("#changeset#", changeset[:7])
	header = header.replace("#contributors#", contributors)
	return header.replace("\t", "").replace("!!", "   ! ")

#################################
# breaks a string (mystr), in piece
# of length (sublen), using a
# separator (sep)
def trunc(mystr,sublen,sep):
	astr = mystr.split(sep)
	s = z = ""
	for x in astr:
		z += x + sep
		s += x + sep
		if len(z) > sublen:
			z=""
			s+="\n"
	return s

#################################
#Returns the implicit ode loop
# arguments are number or reactants (nr)
# an products (np)
def get_implicit_ode(nr=3,np=4):
	s = """
n(idx_dummy) = 1.d0
n(idx_g) = 1.d0
n(idx_CR) = 1.d0
"""

	s += "do i=1,nrea\n"

	#r1=arr_r1(i)
	for i in range(nr):
		s += "r"+str(i+1)+" = arr_r"+str(i+1)+"(i)\n"
	#p1=arr_rp(i)
	for i in range(np):
		s += "p"+str(i+1)+" = arr_p"+str(i+1)+"(i)\n"
	s += "rr = k(i)*"+("*".join(["n(r"+str(i+1)+")" for i in range(nr)]))+"\n"

	#dn(r1)=dn(r1)-rr
	for i in range(nr):
		s += "dn(r"+str(i+1)+") = dn(r"+str(i+1)+") - rr\n"
	#dn(p1)=dn(p1)+rr
	for i in range(np):
		s += "dn(p"+str(i+1)+") = dn(p"+str(i+1)+") + rr\n"

	s += "end do\n"
	return s

###########################
#this function returns if the string line
# starts with one of the items in the array
# of string aarg
def lbeg(aarg,line):
	for arg in aarg:
		if line[:len(arg)] == arg: return True
	return False

###########################
#this function returns if the string line
# ends with one of the items in the array
# of string aarg
def lend(aarg,line):
	for arg in aarg:
		if line[len(line)-len(arg):] == arg: return True
	return False

###############################
#reduce the length of the lines
def cutlines(all_row):
	#check line length
	maxlen = 80
	arowl = []
	for line in all_row:
		if len(line.split("!")[0]) > maxlen:
			erow = line.strip().split(",") #explode using commma
			sall = "" #string
			subrow1 = [] #first line
			subrow2 = [] #second line
			#loop over comma separated parts
			for i in range(len(erow)):
				sall += erow[i]+"," #concatenate string
				#append to the first or the second line depending on the length
				if len(sall)>maxlen-4:
					subrow2.append(erow[i])
				else:
					subrow1.append(erow[i])

			#join first and second line
			if subrow1.join().strip() != "":
				arowl.append((",".join(subrow1)).strip()+", &\n")
			if subrow2.join().strip() != "":
				arowl.append((",".join(subrow2))+"\n")
		else:
			arowl.append(line) #append lines of regular length
	return arowl


###############################
#this function indent f90 file and remove multiple blank lines
def indentF90(filename):
	import os
	#check if the file exists else return
	if not os.path.isfile(filename): return

	#open file for indent
	fh = open(filename)
	arow = [] #array of the lines of the indented file
	is_blank = is_amper = False #flags
	nind = 0 #number indent level
	nspace = 2 #number of space for indent
	tokenclose = ["end do","end if","end function","end subroutine","else if","elseif","else","enddo","end module","endif"]
	tokenclose += ["contains","endfunction","endsubroutine","endmodule","end program", "endprogram"]
	tokenopen = ["do ","function","subroutine","contains","else","else if","elseif","module","program"]
	module_head = "!############### MODULE ##############" #module header comment
	module_head_found = False
	for row in fh:
		srow = row.strip() #trim the row
		if module_head in srow: module_head_found = True #do not duplicate module header
		#check module begin
		if lbeg(["module"], srow) and not module_head_found:
			arow.append("\n") #blank line
			arow.append(module_head+"\n") #comment
		if lbeg(tokenclose, srow): nind -= 1 #check if the line ends with one of tokenclose
		indent = (" "*(nind*nspace)) #compute number of spaces for indent
		if is_amper: indent = (" "*(2*nspace)) + indent #increas indent in case of previous &
		if srow.startswith("#"): indent = "" #no indent for pragmas
		if not(srow=="" and is_blank): arow.append(indent+srow+"\n") #append indented line to array of rows
		is_amper = False #is a line after ampersend flag
		if lend(["&"], srow): is_amper = True #check if the line ends with &
		is_blank = (srow=="") #flag for blank line mode
		if lbeg(tokenopen, srow): nind += 1 #check if the line ends with one of tokenclose
		if lbeg(["if"],srow) and "then" in srow: nind += 1 #check if line stats with if and has then
		if srow == "do": nind += 1
	fh.close()

	arowl = arow[:]

	#write the new file
	fh = open(filename,"w")
	for x in arowl:
		fh.write(x.rstrip()+"\n")
	fh.close()

################################
#This function writes an error
# and exit
def die(msg):
	import sys
	print(msg)
	sys.exit()

#################################
#This function returns a random quotation properly formtatted
# if qall=True print all the quotes
def get_quote(qall=False):
	import random
	quotes = [["If you lie to the computer, it will get you.","Perry Farrar"],
	["Premature optimization is the root of all evil.","Donald Knuth"],
	["Computers are good at following instructions, but not at reading your mind.","Donald Knuth"],
	["Computer Science is embarrassed by the computer.","Alan Perlis"],
	["Prolonged contact with the computer turns mathematicians into clerks and vice versa.","Alan Perlis"],
	["There are two ways to write error-free programs; only the third one works.","Alan Perlis"],
	["Software and cathedrals are much the same - first we build them, then we pray.","Sam Redwine"],
	["Estimate always goes wrong.","Sumit Agrawal"],
	["Weinberg's Second Law: If builders built buildings the way programmers wrote programs, then the first woodpecker that came"\
	 +" along would destroy civilization.","Gerald Weinberg"],
	["Any sufficiently advanced magic is indistinguishable from a rigged demonstration.",""],
	["Any given program, when running, is obsolete.",""],
	["Programming would be so much easier without all the users.",""],
	["Your Zip file is open.",""],
	["Testing can only prove the presence of bugs, not their absence.","Edsger W. Dijkstra"],
	["If debugging is the process of removing bugs, then programming must be the process of putting them in.","Edsger W. Dijkstra"],
	["God is Real, unless declared Integer.","J. Allan Toogood"],
	["Curiously enough, the only thing that went through the mind of the bowl of petunias as it fell was Oh no, not again.","The"\
	 +" Hitchhiker's Guide to the Galaxy"],
	["Computer science differs from physics in that it is not actually a science.","Richard Feynman"],
	["The purpose of computing is insight, not numbers.","Richard Hamming"],
	["Computer science is neither mathematics nor electrical engineering.","Alan Perlis"],
	["I can't be as confident about computer science as I can about biology. Biology easily has 500 years of exciting problems to work"\
	 +" on. It's at that level.","Donald Knuth"],
	["The only legitimate use of a computer is to play games.","Eugene Jarvis"],
	["UNIX is user-friendly, it just chooses its friends.","Andreas Bogk"],
	["Quantum mechanic Seth Lloyd says the universe is one giant, hackable computer. Let's hope it's not running Windows.","Kevin Kelly"],
	["Computers are useless. They can only give you answers.","Pablo Picasso"],
	["Computers in the future may weigh no more than 1.5 tons.","Popular Mechanics (1949)"],
	["Don't trust a computer you can't throw out a window.","Steve Wozniak"],
	["Computers are like bikinis. They save people a lot of guesswork.","Sam Ewing"],
	["If the automobile had followed the same development cycle as the computer, a Rolls-Royce would today cost $100, get a million"\
	 +" miles per gallon, and explode once a year, killing everyone inside.","Robert X. Cringely"],
	["Computers are getting smarter all the time. Scientists tell us that soon they will be able to talk to us.  (And by 'they',"\
	 " I mean 'computers'.  I doubt scientists will ever be able to talk to us.)","Dave Barry"],
	["Most software today is very much like an Egyptian pyramid with millions of bricks piled on top of each other, with no structural\
	 integrity, but just done by brute force and thousands of slaves.","Alan Kay"],
	["No matter how slick the demo is in rehearsal, when you do it in front of a live audience, the probability of a flawless "\
	 + "presentation is inversely proportional to the number of people watching, raised to the power of the amount of money involved.",\
	"Mark Gibbs"],
	["Controlling complexity is the essence of computer programming.","Brian Kernigan"],
	["Software suppliers are trying to make their software packages more 'user-friendly'...  Their best approach so far has been to take"\
	 +" all the old brochures and stamp the words 'user-friendly' on the cover.","Bill Gates"],
	["Programmers are in a race with the Universe to create bigger and better idiot-proof programs, while the Universe is trying to"\
	 + " create bigger and better idiots. So far the Universe is winning.","Rich Cook"],
	["To iterate is human, to recurse divine.","L. Peter Deutsch"],
	["Should array indices start at 0 or 1?  My compromise of 0.5 was rejected without, I thought, proper consideration.","Stan Kelly-Bootle"],
	["Any code of your own that you haven't looked at for six or more months might as well have been written by someone else.","Eagleson's Law"],
	["All science is either physics or stamp collecting.", "Ernest Rutherford"],
	["Done is better than perfect.", ""],
	["Computers are like Old Testament gods; lots of rules and no mercy","Joseph Campbell"],
	["A computer lets you make more mistakes faster than any other invention with the possible exceptions of handguns and Tequila.",\
	 "Mitch Ratcliffe"],
	["Computer Science is no more about computers than astronomy is about telescopes.","Edsger W. Dijkstra"],
	["To err is human, but to really foul things up you need a computer.","Paul Ehrlich"],
	["Debugging is twice as hard as writing the code in the first place. Therefore, if you write the code as cleverly as possible, you are,"\
	 +" by definition, not smart enough to debug it.","Brian W. Kernighan"],
	["Always code as if the guy who ends up maintaining your code will be a violent psychopath who knows where you live.","Martin Golding"],
	["One of my most productive days was throwing away 1000 lines of code.","Ken Thompson "],
	["And God said, \"Let there be light\" and segmentation fault (core dumped)",""],
	["Today, most software exists, not to solve a problem, but to interface with other software","Ian Angell"],
	["Measure twice, cut once",""],
	["Weeks of programming can save you hours of planning",""],
	["All models are wrong; some models are useful","George Box"],
	["The generation of random numbers is too important to be left to chance","Robert Coveyou"],
	["Problems worthy / of attack / prove their worth / by hitting back","Piet Hein"],
	["Good, Fast, Cheap: Pick any two","Memorandum RFC 1925"],
	["One size never fits all","Memorandum RFC 1925"],
	["No matter how hard you push and no matter what the priority,you can't increase the speed of light","Memorandum RFC 1925"],
	["Chemistry has been termed by the physicist as the messy part of physics", "Frederick Soddy "],
	["Don't worry if it doesn't work right. If everything did, you'd be out of a job.", "Mosher's Law"],
	["Beware of bugs in the above code; I have only proved it correct, not tried it.", "Donald Knuth"],
	["Given enough eyeballs, all bugs are shallow.", "Eric S. Raymond"],
	["We make an extreme, but wholly defensible, statement: There are no good, general methods for solving systems of more than one nonlinear equation.",\
	 "Numerical Recipes in C"]
	]
	qrange = 1
	print("")
	if(qall): qrange = len(quotes)
	for i in range(qrange):
		irand = int(random.random()*(len(quotes)))
		if(qall): irand = i
		qtup = quotes[irand]
		myqt = trunc(str(irand+1)+". "+qtup[0],40," ").upper().strip()
		amyqt = myqt.split("\n")
		lqt = max([len(x) for x in amyqt])
		print("")
		if i == 0: print("*"*lqt)
		print(myqt)
		if qtup[1].strip() == "": qtup[1] = "Anonymous"
		print("--- "+qtup[1])
		if i == qrange - 1:
			print("*"*lqt)
			print("")


def keyb_input(message):
	try:
		raw_input(message)
	except:
		input(message)
