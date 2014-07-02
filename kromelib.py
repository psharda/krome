#KROME is a nice and friendly chemistry package for a wide range of 
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
# Written and developed by Tommaso Grassi
# tommasograssi@gmail.com,
# Starplan Center, Copenhagen.
# Niels Bohr Institute, Copenhagen.
#
# Co-developer Stefano Bovino
# sbovino@astro.physik.uni-goettingen.de
# Institut fuer Astrophysik, Goettingen.
#
# Others (alphabetically): F.A. Gianturco, J.Prieto,
# D.R.G. Schleicher, D. Seifried, E. Simoncini 
#
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
##################################
class molec():
	name = "" #molecule name
	charge = 0 #charge (0=neutral)
	mass = 0. #mass in g
	zatom = 0 #atomic number (also for molecules)
	neutrons = 0 #total number of neutrons
	ename = "" #exploded name (e.g. H2C4=CCCCHH)
	fname = "" #f90 name (e.g. H+=Hj, D-=Dk)
	phname = "" #compact name for photoionizations (e.g. Fe++++=Fej4)
	fidx = "idx_" #index for fortran (e.g. H+=idx_Hj, D-=idx_Dk)
	is_atom = False #flag to identify atoms
	chempot = 0. #chemical potential (J/mol)
	poly1 = [0.e0]*7 #nasa polynomials (200-1000K)
	poly2 = [0.e0]*7 #nasa polynomials (1000-5000K)
	Tpoly = [0.e0]*3 #temperature limits
	idx = 0 #species index
	enthalpy = 0.e0 #enthalpy of formation
	atomcount = dict() #dictionary containin the count of atoms including zero (e.g H2O is {"H":2, "O":1})
	atomcount2 = dict() #dictionary containin the count of atoms without zero species (e.g H2O is {"H":2, "O":1})
	natoms = 0 #the number of atoms (e.g. diatomic=2)
	ve_vib = "__NONE__" #vibrational constant in K
	be_rot = "__NONE__" #rotational constant in K

	def __init__(self):
		self.poly1 = [0.e0]*7
		self.poly2 = [0.e0]*7
		self.Tpoly = [0.e0]*3
		self.atomcount = dict()	
##################################
class reaction():
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
	ifrate = "" #if condition on rate, e.g. if(Tgas>1d2):
	isCR = False #flag this reaction as CR
	isXRay = False #flag this reaction as XRay
	#method: constructor to initialize lists
	def __init__(self):
		self.reactants = []
		self.products = []

	#method: build verbatim from reactants and products
	def build_verbatim(self):
		myr = []
		myp = []
		for r in self.reactants:
			if(r.name!="dummy"):
				myr.append(r.name)	
		for p in self.products:
			if(p.name!="dummy"):
				myp.append(p.name)
		self.verbatim = " + ".join(myr)+" -> "+" + ".join(myp)

	#method: build photochemical rate
	def build_phrate(self,photoBlock=False):
		if(not("krome_kph_auto" in self.krate) and not(photoBlock)): return
		myr = self.reactants
		self.kphrate = self.krate.replace("krome_kph_auto=","")
		self.krate = "photoBinRates("+str(self.idxph)+")"
		#if(self.krate.strip()=="krome_kph_auto"): 
		#	self.krate = "krome_kph_"+myr[0].phname
		#else: 
		#	self.krate = "krome_kph_" + myr[0].phname.capitalize() + "R" + str(self.idx)

	#method: build RHS
	def build_RHS(self,useNuclearMult=False):
		if(self.idx<=0):
			print "************************************************************"
			print "ERROR: reaction index (reaction.idx) must be greater than 0"
			print "Probably you have to define idx"
			print "************************************************************"
			sys.exit()

		#keeps into account nuclear multeplicity (if option is enabled)
		nuclearMult = "" #default nuclear multeplicity value
		if(useNuclearMult):
			uncurledR = [self.reactants[i].name for i in range(len(self.reactants)) if not(self.curlyR[i])]
			if(len(uncurledR)==2):
				if(uncurledR[0]==uncurledR[1]): nuclearMult = "0.5d0*"
			if(len(uncurledR)==3):
				if(uncurledR[0]==uncurledR[1] and uncurledR[1]==uncurledR[2]): nuclearMult = "0.16666666667d0*"
				if(uncurledR[0]==uncurledR[1] and uncurledR[1]!=uncurledR[2]): nuclearMult = "0.5d0*"
				if(uncurledR[1]==uncurledR[2] and uncurledR[0]!=uncurledR[1]): nuclearMult = "0.5d0*"
		self.nuclearMult = nuclearMult
		self.RHS = nuclearMult+"k("+str(self.idx)+")"
		self.RHSvar = "kflux"+str(self.idx)
		ns = []
		i = 0
		for r in self.reactants:
			if(self.curlyR[i]): continue #skip curly reactants
			i += 1
			ns.append("n("+str(r.fidx)+")")
			if(r.idx<=0):
				print "************************************************************"
				print "ERROR: species index (molec.idx) must be greater than 0"
				print "Probably you have to define idx"
				print "************************************************************"
				sys.exit()
		if(len(ns)>0): self.RHS += "*" + ("*".join(ns))

	#method:build pseudo hash (for unique reactions to avoid duplicates)
	def build_pseudo_hash(self):
		rname = []
		pname = []
		for r in self.reactants:
			rname.append(r.name)
		
		for p in self.products:
			pname.append(p.name)
		self.pseudo_hash = ("_".join(sorted(rname)))+"|"+("_".join(sorted(pname))) 

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
		if(mass_reactants!=0 and (mode=="ALL" or "MASS" in mode)):
			if(abs(1.e0-mass_products/mass_reactants)>1e-6):
				print "************************************************"
				print "WARNING: problem with mass conservation in reaction", self.idx
				print "reaction:",self.verbatim
				print "reactants:", [r.name for r in self.reactants], mass_reactants
				print "products:", [p.name for p in self.products], mass_products
				print "mass ratio (prods/reacts):",mass_products/mass_reactants, "(should be 1.0)"
				print "You can remove this check with the -nomassCheck option"	
				print "************************************************"
				a = raw_input("Any key to continue q to quit... ")
				if(a=="q"): print sys.exit()
		if(abs(charge_products - charge_reactants)!=0 and (mode=="ALL" or "CHARGE" in mode)):
			print "************************************************"
			print "WARNING: problem with charge conservation in reaction", self.idx
			print "reaction:",self.verbatim
			print "reactants:", [r.name for r in self.reactants], charge_reactants
			print "products:", [p.name for p in self.products], charge_products
			print "You can remove this check with the -nochargeCheck option"
			print "************************************************"
			a = raw_input("Any key to continue q to quit... ")
			if(a=="q"): print sys.exit()

	#calcluate enthalpy of formation
	def enthalpy(self):
		if("krome_kph" in self.krate): return

		if(len(self.reactants)==0 and len(self.products)==0):
			self.dH = 0e0
			return
		if(len(self.reactants)==0 or len(self.products)==0):
			print "ERROR: you have called enthalpy calculation"
			print " with empty reactants and/or products!"
			sys.exit()
		reag = self.reactants #copy reactants
		prod = self.products #copy products
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
		if(available): self.dH = (rH-pH)*1.60217657e-12 #eV->erg (cooling<0)


	#calculate reverse reaction using polynomials
	def doReverse(self):
		pidx = "(/"+(",".join([x.fidx for x in self.products]))+"/)"
		ridx = "(/"+(",".join([x.fidx for x in self.reactants]))+"/)"
		ndif = len(self.reactants)-len(self.products)
		kk = "("+self.krate+") / exp(revKc(Tgas,"+ridx+","+pidx+"))"
		if(ndif!=0): kk ="0.d0"  #" * (1.3806488d-2 * Tgas)**("+str(ndif)+")"
		return kk

#################################
#create tabvar (probably not the best interface ever)
def create_tabvar(mytabvar,mytabpath,mytabxxyy,anytabvars,anytabfiles,anytabpaths,anytabsizes,coevars,ivarcoe):
	if(mytabvar.split("_")[0].lower()!="user"):
		print "ERROR: to avoid conflicts common variables with @tabvar should begin with user_"
		print " you provided: "+mytabvar
		print " it should be: user_"+mytabvar
		sys.exit()
	#check if file exists
	if(not(file_exists(mytabpath))):
		print "ERROR: file "+mytabpath+" not found!"
		print " note that the path must be relative to the ./krome command"
		sys.exit()

	#read the size of the table from the first line file
	fhtab = open(mytabpath,"rb")
	for tabrow in fhtab:
		stabrow = tabrow.strip()
		if("," in stabrow): 
			mytabsize = [xx.strip() for xx in stabrow.split(",")]
		else:
			print "ERROR: the file "+mytabpath+" must contain the size of the"
			print " table in the first line (comma separated, e.g. 50,30)"
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
	tabf =  "fit_anytab2D("+anytabx+", &\n"+anytaby+", &\n"+anytabz+", &\n"+anytabxmul+", &\n"+anytabymul+", &\n"+mytabxxyy+")"
	coevars[mytabvar] = [ivarcoe,tabf]

	print "Found tabvar:",mytabvar,"("+mytabpath+")", "["+(",".join(mytabsize))+"]"

#############################
def addVarCoe(mytabvar,tabf,coevars,ivarcoe):
	coevars[mytabvar] = [ivarcoe,tabf]

#############################
#cooling index list
def get_cooling_index_list():
	idxcoo = {"H2":1,"H2GP":2,"atomic":3, "CEN":3, "HD":4, "Z":5, "metal":5, "dH":6, "enthalpic":6, "dust":7,\
		"compton":8,"CIE":9, "continuum":10, "cont":10,"exp":11,"expansion":11}

	idxscoo = []
	maxv = 0
	for (k,v) in idxcoo.iteritems():
		idxscoo.append([v,"idx_cool_"+k+" = "+str(v)])
		maxv = max(maxv,v)
	idxscoo = sorted(idxscoo,key=lambda x:x[0])
	idxscoo.append([99,"ncools = "+str(maxv)])
	return [x[1] for x in idxscoo]

#############################
#heating index list
def get_heating_index_list():
	idxhea = {"chem":1,"compress":2, "compr":2, "photo":3, "dH":4, "enthalpic":4, "photoAv":5, "Av":5,\
		"CR":6, "dust":7, "xray":8}

	idxshea = []
	maxv = 0
	for (k,v) in idxhea.iteritems():
		idxshea.append([v, "idx_heat_"+k+" = "+str(v)])
		maxv = max(maxv,v)
	idxshea = sorted(idxshea,key=lambda x:x[0])
	idxshea.append([99,"nheats = "+str(maxv)])
	return [x[1] for x in idxshea]


####################################
#solar metallicities
def get_solar_abundances():
	#solar abundances from Anders+Grevesse 1989
	solar_abs = {
		"Li":"2.046595d-9",
		"C" :"3.620072d-4",
		"N" :"1.121864d-4",
		"O" :"8.530466d-4",
		"F" :"3.021505d-8",
		"Ne":"1.232975d-4",
		"Na":"2.057348d-6",
		"Mg":"3.849462d-5",
		"Al":"3.043011d-6",
		"Si":"3.584229d-5",
		"P" :"3.727599d-7",
		"S" :"1.845878d-5",
		"Cl":"1.878136d-7",
		"Ca":"2.189964d-6",
		"Fe":"3.225806d-5"
		}
	return solar_abs



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
	if(arg in ve):
		return ve[arg]*1.42879e0 #cm-1 to K 
	else:
		False

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
	if(arg in ve):
		return ve[arg]*1.42879e0 #cm-1 to K 
	else:
		False

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
		print "ERROR: expected integer, got " + type(input)
		sys.exit()
	if not 0 < input < 4000:
		print "ERROR: Argument must be between 1 and 3999"
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
def get_example(nsp,useX):
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
	if(useX):
		sfile = sfile.replace("@rho@",",rho").replace("@nH@","1.d0").replace("@norm@","x(:) = x(:) / sum(x) !normalize")
		sfile = sfile.replace("@rho_init@","rho = 1d-18 !gas density (g/cm3)").replace("@rhof@"," rho,")
	else:
		sfile = sfile.replace("@rho@","").replace("@nH@","1.d4").replace("@norm@","")
		sfile = sfile.replace("@rho_init@","").replace("@rhof@","")
	sfile = sfile.replace("@nsp@",str(nsp))

	if("@" in sfile):
		print sfile
		print "ERROR: missing replacement in get_example() function!"
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
def get_Tshortcut(rea,slist,cvars=[]):
	shcut = ["t = Tgas !alias for Tgas (K)",
	"logT = log10(Tgas) !log10 of Tgas (#)",
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
	krea = rea.krate
	for x in shcut:
		ax = x.split("=") #split the shortcut
		xvar = parsevar(krea) #parse the variable in the rate coefficient
		krea = krea.replace(ax[0].strip(),"")
		if((ax[0].strip().lower() in xvar) and not(x in slist)):
			#search for dependencies between shortcuts
			xtmp = x
			for xx in shcut:
				axx = xx.split("=") #split the shortcut
				xxvar = parsevar(xtmp)  #parse the variables in the rate coeffcient found
				xtmp = xtmp.replace(axx[0].strip(),"")
				if((axx[0].strip().lower() in xxvar) and not(xx in slist)):
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
			if(xkey.lower().strip()==cv.lower().strip()): 
				xFound = True
				break
		if(not(xFound)):
			slistu.append(x)

	#keep only unique shortcuts (remove duplicates)
	ulist = []
	for x in slistu:
		if(not(x in ulist)): ulist.append(x)

	return ulist

##################################
#get list of available commands
def get_usage():
	print "use -h to see the help"
	sys.exit()

##################################
#truncate F90 expression using sep as separator for blocks shorter than sublen
def truncF90(mystr, sublen, sep):
	#split (&\n) the string mystr in parts smaller tha sublen using sep as separator
	if(mystr.strip()==""): return mystr
	mystr = mystr.replace("**","##")
	astr = mystr.split(sep)
	s = z = ""
	first = True
	for x in astr:
		if(len(z+x)>sublen and not(first)):
			s += "&\n"
			z = ""
		zep = sep
		if(first): zep = "" 
		s += zep + x
		z += zep + x
		first = False
	return s.replace("##","**")	

################################
def at_extract(arg):
	aarg = arg.split(":")
	if(len(aarg)!=2):
		print "ERROR: @label:value format not respected for"
		print arg
		sys.exit()
	aarg[0] = aarg[0].replace("@","")
	return {aarg[0]:aarg[1]}

##################################
#double format for f90 expressions
def format_double(snum):
	snum = str(snum)
	#format string number to F90 double
	if("d" in snum): return snum
	if("e" in snum): return snum.replace("e","d")	
	return snum+"d0"

##################################
#format the subelement according to its size (Fe, N, ...)
def format_subel(subel):	
	if(len(subel)==1):
		fsubel = subel.upper()
		if(fsubel=="G"): fsubel = "g"
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
def parser(name, mass_dic, atoms, thermo_data):

	mymol = molec() #oggetto molec
	namecp = name.upper()
	if(namecp=="E-"): namecp = "E" #avoid double negative charge
	ename = [] #exploded name
	mass = 0. #init mass
	is_atom = True #atom flag
	founds = 0 #atoms found

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
		"TI":22,
		"CA":20,
		"CR":24,
		"MN":25,
		"FE":26,
		"NI":28}


	#add isotopes Z numbers (same as non-isotopes)
	zdic_copy = zdic.copy()
	for k,v in zdic_copy.iteritems():
		for i in range(2,56):
			zdic[str(i)+k] = v

	#check for fake species
	if("FK" in name):
		mymol.name = name #name
		mymol.mass = 0e0 #mass (g)
		mymol.ename = name #exploded name
		mymol.charge = 0 #charge
		mymol.zatom = 0 #atomic number
		mymol.fname = name #f90 name
		mymol.is_atom = True #atom flag
		mymol.fidx = "idx_"+name #f90 index
		mymol.neutrons = 0 #number of neutrons
		return mymol
	
	zatom = 0 #atomic number init
	#loop over charcters
	for atm in atoms:
		a = atm.upper() #capitalize name
		if(not(a) in namecp): continue #skip
		#loop to up to _30 subscript
		for j in range(30):
			if(a in namecp):
				idx = namecp.find(a) #find position
				subs = a
				mult = "0" #multiplicator
				for i in range(idx+len(a),len(namecp)):
					if(not(is_number(namecp[i]))): break
					mult += namecp[i] #find multiplicator
					subs += namecp[i]
				imult = max(int(mult),1) #evaluate multiplicator (must be >0)
				mass += mass_dic[a]*imult #compute mass
				if(a in zdic): zatom += zdic[a]*imult #increase atomic number
				if(format_subel(a) in mymol.atomcount): 
					mymol.atomcount[format_subel(a)] += imult #increase atom count
				else:
					mymol.atomcount[format_subel(a)] = imult #init atom count
				ename += [format_subel(a)]*imult #exploded name
				namecp = namecp.replace(subs,"",1) #remove found in name
				if(namecp==""): break #if nothing more to find break loop
		if(a!="+" and a!="-"): founds += imult #count found atoms for is_atom
	if(founds>1): is_atom = False #atoms have only one atom (viz.)
	
	mymol.atomcount2 = dict()
	natoms = 0
	for (k,v) in mymol.atomcount.iteritems():
		if(v>0):
			mymol.atomcount2[k] = v
			if(k!="+" and k!="-"): natoms += v

	mymol.natoms = natoms #number of atoms (e.g. diatom=2)

	#get vibrational constant in K
	if(get_ve_vib(name)):
		mymol.ve_vib = get_ve_vib(name)
	#get rotational constant in K
	if(get_be_rot(name)):
		mymol.be_rot = get_be_rot(name)

	mymol.name = name #name
	mymol.mass = mass #mass (g)
	mymol.ename = sorted(ename) #exploded name
	mymol.charge = 0 #charge
	mymol.zatom = zatom #atomic number
	mymol.fname = name.replace("+","j").replace("-","k") #f90 name
	mymol.is_atom = is_atom #atom flag
	f90idx = "idx_"+name.replace("+","j").replace("-","k").replace("(","_").replace(")","").replace("[","").replace("]","_") #f90 index
	if(f90idx.endswith("_")): f90idx = f90idx[:-1] #remove last underscore if any
	mymol.fidx = f90idx #index in f90 format
	
	if("+" in name): mymol.charge = name.count("+") #get + charge
	if("-" in name): mymol.charge = -name.count("-") #get - charge

	#number of neutrons (computed using total mass)
	Nn = round((mymol.mass - (me*(mymol.zatom - mymol.charge) + mp*(mymol.zatom))) / mn,0)
	mymol.neutrons = int(Nn)


	#name for photoionization reactions (e.g. Sijjj = Sij3)
	jj = kk = ""
	if(mymol.charge==1): jj = "j"
	if(mymol.charge>1): jj = "j"+str(mymol.charge)
	if(mymol.charge==-1): kk = "k"
	if(mymol.charge<-1): kk = "k"+str(mymol.charge)
	mymol.phname = name.replace("+","").replace("-","") + jj + kk

	#electron has negative charge
	if(mymol.name=="E"): mymol.charge = -1
	
	#thermal data
	if(mymol.name in thermo_data):
		mymol.poly1 = thermo_data[mymol.name][10:] #NASA polynomials lower T interval (min-med)
		mymol.poly2 = thermo_data[mymol.name][3:10] #NASA polynomials upper T interval (med-max)
		mymol.Tpoly = thermo_data[mymol.name][0:3] #(K) [min,med,max] T interval limits

	#compute enthaly @300K using NASA poly
	p = mymol.poly1 #copy polynomials
	Tgas = 300. #K
	polyH = p[0] + p[1]*0.5*Tgas + p[2]*Tgas**2/3. + p[3]*Tgas**3*0.25 + p[4]*Tgas**4*0.2 + p[5]/Tgas
	mymol.enthalpy = polyH*8.314472e-3*Tgas*0.01036410e0 #eV

	#checks parsing results
	if(len(namecp)>0): 
		print "************************************************"
		print "ERROR: Parsing problem for", name
		print "Unknown subelements in substring \"" + namecp +"\"."
		print "Probably you have to add some subelements to the dictionary mass_dic."
		if(len(atoms)<30):
			print "Dictionary now contains the following subelements:"
			print atoms
		print "************************************************"
		sys.exit()
	return mymol

######################################
def get_file_list():
	files = []
	files.append("kromelib.py")
	files.append("patches")
	files.append("patches/flash")
	files.append("patches/flash/Driver")
	files.append("patches/flash/Driver/DriverMain")
	files.append("patches/flash/Driver/DriverMain/Driver_sourceTerms.F90")
	files.append("patches/flash/Driver/DriverMain/Driver_initSourceTerms.F90")
	files.append("patches/flash/Driver/DriverMain/Driver_finalizeSourceTerms.F90")
	files.append("patches/flash/physics")
	files.append("patches/flash/physics/sourceTerms")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/Makefile")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistry_finalize.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistry_init.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistry.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/pchem_mapNetworkToSpecies.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/Makefile")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/Config")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/KromeChemistry_init.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/KromeChemistry.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistryMain/KromeChemistry_data.F90")
	files.append("patches/flash/physics/sourceTerms/KromeChemistry/KromeChemistry_interface.F90")
	files.append("patches/flash/Simulation")
	files.append("patches/flash/Simulation/SimulationMain")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/flash.par")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/Simulation_data.F90")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/Makefile")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/Simulation_init.F90")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/flash.par_1d")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/Config")
	files.append("patches/flash/Simulation/SimulationMain/Chemistry_Krome_Collapse/Simulation_initBlock.F90")
	files.append("patches/flash/Simulation/SimulationComposition")
	files.append("patches/flash/Simulation/SimulationComposition/KromeChemistry")
	files.append("patches/flash/Simulation/SimulationComposition/KromeChemistry/SpeciesList.txt")
	files.append("patches/flash/Simulation/SimulationComposition/KromeChemistry/Config")
	files.append("patches/flash/Simulation/SimulationComposition/KromeChemistry/Simulation_initSpecies.F90")
	files.append("patches/enzo")
	files.append("patches/enzo/kromebuild.sh")
	files.append("patches/enzo/InitializeRateData.C")
	files.append("patches/enzo/Grid_SolveRateAndCoolEquations.C")
	files.append("patches/enzo/krome_initab.F90")
	files.append("patches/enzo/evaluate_tgas.F90")
	files.append("patches/enzo/Make.config.objects")
	files.append("patches/enzo/Grid_IdentifySpeciesFieldsKrome.C")
	files.append("patches/enzo/krome_driver.F90")
	files.append("patches/enzo/Grid.h")
	files.append("patches/ramses")
	files.append("patches/ramses/condinit.f90")
	files.append("patches/ramses/init_flow_fine.f90")
	files.append("patches/ramses/output_hydro.f90")
	files.append("patches/ramses/hydro_parameters.f90")
	files.append("patches/ramses/amr_parameters.f90")
	files.append("patches/ramses/Makefile")
	files.append("patches/ramses/cooling_fine.f90")
	files.append("patches/ramses/cooling_module.f90")
	files.append("patches/ramses/read_hydro_params.f90")
	files.append("tests")
	files.append("tests/test.f90")
	files.append("tests/Makefile")
	files.append("tests/MakefileF90")
	files.append("tests/compact")
	files.append("tests/compact/test.f90")
	files.append("tests/compact/shock1D.dat")
	files.append("tests/compact/MakefileCompact")
	files.append("tests/compact/plot.gps")
	files.append("tests/map")
	files.append("tests/map/test.f90")
	files.append("tests/map/plot.gps")
	files.append("tests/collapse")
	files.append("tests/collapse/test.f90")
	files.append("tests/collapse/plot.gps")
	files.append("tests/cloud")
	files.append("tests/cloud/test.f90")
	files.append("tests/cloud/Makefile")
	files.append("tests/cloud/plot.gps")
	files.append("tests/collapseUV")
	files.append("tests/collapseUV/test.f90")
	files.append("tests/collapseUV/plot.gps")
	files.append("tests/collapseZ")
	files.append("tests/collapseZ/test.f90")
	files.append("tests/collapseZ/MakefileF90")
	files.append("tests/collapseZ/plot.gps")
	files.append("tests/slowmanifold")
	files.append("tests/slowmanifold/test.f90")
	files.append("tests/slowmanifold/testzzz.f90")
	files.append("tests/slowmanifold/plot.gps")
	files.append("tests/atmosphere")
	files.append("tests/atmosphere/test.f90")
	files.append("tests/atmosphere/eddy.dat")
	files.append("tests/atmosphere/init_spec.dat")
	files.append("tests/atmosphere/layers_data")
	files.append("tests/atmosphere/plot.gps")
	files.append("tests/MakefileCompact")
	files.append("tests/shock1Dcool")
	files.append("tests/shock1Dcool/test.f90")
	files.append("tests/shock1Dcool/shock1D.dat")
	files.append("tests/shock1Dcool/plot.gps")
	files.append("tests/dust")
	files.append("tests/dust/test.f90")
	files.append("tests/dust/plot.gps")
	files.append("tests/reverse")
	files.append("tests/reverse/test.f90")
	files.append("tests/reverse/plot.gps")
	files.append("tests/shock1D")
	files.append("tests/shock1D/test.f90")
	files.append("tests/shock1D/plot.gps")
	files.append("tests/lotkav")
	files.append("tests/lotkav/test.f90")
	files.append("tests/lotkav/plot.gps")
	files.append("tests/lotkav/lotkav")
	files.append("tests/Makefile_pedantic")
	files.append("tests/wrapC")
	files.append("tests/wrapC/test.c")
	files.append("tests/wrapC/krome.h")
	files.append("tests/shock1Dphoto")
	files.append("tests/shock1Dphoto/test.f90")
	files.append("tests/shock1Dphoto/shock1D.dat")
	files.append("tests/shock1Dphoto/plot.gps")
	files.append("krome")
	files.append("kromeobj.py")
	files.append("build")
	files.append("tools")
	files.append("src")
	files.append("src/krome_constants.f90")
	files.append("src/krome_tabs.f90")
	files.append("src/kromeF90.f90")
	files.append("src/krome_user.f90")
	files.append("src/krome_stars.f90")
	files.append("src/krome_user_commons.f90")
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
	files.append("data/thermo30.dat")
	files.append("data/coolZ.dat")
	files.append("solver")
	files.append("solver/nleq_all.f")
	files.append("solver/opkda2.f")
	files.append("solver/opkdmain.f")
	files.append("solver/dvode_f90_license.txt")
	files.append("solver/opkda1.f")
	files.append("solver/dvode_f90_m.f90")
	return files
	files.append("networks")
	files.append("networks/react_primordial2")
	files.append("networks/react_primordial_photo")
	files.append("networks/react_NO")
	files.append("networks/react_kast80")
	files.append("networks/react_SM")
	files.append("networks/react_primordial_photoH2")
	files.append("networks/react_primordialZ")
	files.append("networks/react_lowmetal")
	files.append("networks/react_primordialZ2")
	files.append("networks/react_primordial")
	files.append("networks/react_cloud")
	files.append("networks/react_star")
	files.append("networks/react_dummy")
	files.append("networks/react_primordial_UV")

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
	    except Exception, e:
		print e

#################################
#verner96 photochemistry cross sections dictionary
def get_photo_crossV96(atom):
	myatom = atom.lower()
	cross = {"h" : "sigma_h = sigma_v96(energy_eV,4.298d-01,5.475d+04,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"he" : "sigma_he = sigma_v96(energy_eV,1.361d+01,9.492d+02,1.469d+00,3.188d+00,2.039d+00,4.434d-01,2.136d+00)",
		"hej" : "sigma_hej = sigma_v96(energy_eV,1.720d+00,1.369d+04,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"li" : "sigma_li = sigma_v96(energy_eV,3.107d+00,6.245d+01,1.501d+01,4.895d+00,0.000d+00,0.000d+00,0.000d+00)",
		"lij" : "sigma_lij = sigma_v96(energy_eV,2.006d+01,3.201d+02,7.391d+00,2.916d+00,0.000d+00,0.000d+00,0.000d+00)",
		"lij2" : "sigma_lij2 = sigma_v96(energy_eV,3.871d+00,6.083d+03,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"be" : "sigma_be = sigma_v96(energy_eV,9.539d+00,2.932d+05,4.301d-01,1.052d+01,3.655d-01,8.278d-04,1.269d-02)",
		"bej" : "sigma_bej = sigma_v96(energy_eV,1.181d+00,2.678d+02,5.645d+00,1.170d+01,0.000d+00,0.000d+00,0.000d+00)",
		"bej2" : "sigma_bej2 = sigma_v96(energy_eV,1.760d+01,5.458d+02,1.719d+01,3.157d+00,0.000d+00,0.000d+00,0.000d+00)",
		"bej3" : "sigma_bej3 = sigma_v96(energy_eV,6.879d+00,3.422d+03,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"b" : "sigma_b = sigma_v96(energy_eV,5.213d-01,5.466d+00,8.618d+00,1.728d+01,1.887d+01,1.319d+01,4.556d+00)",
		"bj" : "sigma_bj = sigma_v96(energy_eV,2.869d+00,1.859d+04,1.783d+00,1.618d+01,3.503d+00,4.960d-03,3.400d-02)",
		"bj2" : "sigma_bj2 = sigma_v96(energy_eV,1.041d+00,5.393d+01,1.767d+01,9.540d+00,0.000d+00,0.000d+00,0.000d+00)",
		"bj3" : "sigma_bj3 = sigma_v96(energy_eV,3.336d+01,2.846d+02,2.163d+01,2.624d+00,0.000d+00,0.000d+00,0.000d+00)",
		"bj4" : "sigma_bj4 = sigma_v96(energy_eV,1.075d+01,2.190d+03,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"c" : "sigma_c = sigma_v96(energy_eV,2.144d+00,5.027d+02,6.216d+01,5.101d+00,9.157d-02,1.133d+00,1.607d+00)",
		"cj" : "sigma_cj = sigma_v96(energy_eV,4.058d-01,8.709d+00,1.261d+02,8.578d+00,2.093d+00,4.929d+01,3.234d+00)",
		"cj2" : "sigma_cj2 = sigma_v96(energy_eV,4.614d+00,1.539d+04,1.737d+00,1.593d+01,5.922d+00,4.378d-03,2.528d-02)",
		"cj3" : "sigma_cj3 = sigma_v96(energy_eV,3.506d+00,1.068d+02,1.436d+01,7.457d+00,0.000d+00,0.000d+00,0.000d+00)",
		"cj4" : "sigma_cj4 = sigma_v96(energy_eV,4.624d+01,2.344d+02,2.183d+01,2.581d+00,0.000d+00,0.000d+00,0.000d+00)",
		"cj5" : "sigma_cj5 = sigma_v96(energy_eV,1.548d+01,1.521d+03,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"n" : "sigma_n = sigma_v96(energy_eV,4.034d+00,8.235d+02,8.033d+01,3.928d+00,9.097d-02,8.598d-01,2.325d+00)",
		"nj" : "sigma_nj = sigma_v96(energy_eV,6.128d-02,1.944d+00,8.163d+02,8.773d+00,1.043d+01,4.280d+02,2.030d+01)",
		"nj2" : "sigma_nj2 = sigma_v96(energy_eV,2.420d-01,9.375d-01,2.788d+02,9.156d+00,1.850d+00,1.877d+02,3.999d+00)",
		"nj3" : "sigma_nj3 = sigma_v96(energy_eV,5.494d+00,1.690d+04,1.714d+00,1.706d+01,7.904d+00,6.415d-03,1.937d-02)",
		"nj4" : "sigma_nj4 = sigma_v96(energy_eV,4.471d+00,8.376d+01,3.297d+01,6.003d+00,0.000d+00,0.000d+00,0.000d+00)",
		"nj5" : "sigma_nj5 = sigma_v96(energy_eV,6.943d+01,1.519d+02,2.627d+01,2.315d+00,0.000d+00,0.000d+00,0.000d+00)",
		"nj6" : "sigma_nj6 = sigma_v96(energy_eV,2.108d+01,1.117d+03,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"o" : "sigma_o = sigma_v96(energy_eV,1.240d+00,1.745d+03,3.784d+00,1.764d+01,7.589d-02,8.698d+00,1.271d-01)",
		"oj" : "sigma_oj = sigma_v96(energy_eV,1.386d+00,5.967d+01,3.175d+01,8.943d+00,1.934d-02,2.131d+01,1.503d-02)",
		"oj2" : "sigma_oj2 = sigma_v96(energy_eV,1.723d-01,6.753d+02,3.852d+02,6.822d+00,1.191d-01,3.839d-03,4.569d-01)",
		"oj3" : "sigma_oj3 = sigma_v96(energy_eV,2.044d-01,8.659d-01,4.931d+02,8.785d+00,3.143d+00,3.328d+02,4.285d+01)",
		"oj4" : "sigma_oj4 = sigma_v96(energy_eV,2.854d+00,1.642d+04,1.792d+00,2.647d+01,2.836d+01,3.036d-02,5.554d-02)",
		"oj5" : "sigma_oj5 = sigma_v96(energy_eV,7.824d+00,6.864d+01,3.210d+01,5.495d+00,0.000d+00,0.000d+00,0.000d+00)",
		"oj6" : "sigma_oj6 = sigma_v96(energy_eV,8.709d+01,1.329d+02,2.535d+01,2.336d+00,0.000d+00,0.000d+00,0.000d+00)",
		"oj7" : "sigma_oj7 = sigma_v96(energy_eV,2.754d+01,8.554d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"f" : "sigma_f = sigma_v96(energy_eV,1.297d+01,3.803d+03,2.587d+00,7.275d+00,2.170d-03,1.701d-04,1.345d-02)",
		"fj" : "sigma_fj = sigma_v96(energy_eV,1.763d+00,8.013d+01,1.667d+01,1.050d+01,5.103d-01,1.715d+01,7.724d-01)",
		"fj2" : "sigma_fj2 = sigma_v96(energy_eV,2.542d+00,1.541d+02,5.742d+01,6.614d+00,1.115d+00,1.641d+01,5.124d+00)",
		"fj3" : "sigma_fj3 = sigma_v96(energy_eV,7.744d-01,3.165d+00,1.099d+02,9.203d+00,6.812d+00,9.531d+01,9.781d+00)",
		"fj4" : "sigma_fj4 = sigma_v96(energy_eV,7.286d-01,4.690d-01,1.400d+02,9.718d+00,2.570d-01,1.506d+02,2.574d-01)",
		"fj5" : "sigma_fj5 = sigma_v96(energy_eV,4.008d+00,1.157d+04,1.848d+00,2.446d+01,2.411d+01,2.071d-02,3.998d-02)",
		"fj6" : "sigma_fj6 = sigma_v96(energy_eV,2.563d+00,6.930d+01,7.547d+01,6.448d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fj7" : "sigma_fj7 = sigma_v96(energy_eV,1.131d+02,1.039d+02,2.657d+01,2.255d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fj8" : "sigma_fj8 = sigma_v96(energy_eV,3.485d+01,6.759d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"ne" : "sigma_ne = sigma_v96(energy_eV,4.870d+00,4.287d+03,5.798d+00,8.355d+00,2.434d-01,4.236d-02,5.873d+00)",
		"nej" : "sigma_nej = sigma_v96(energy_eV,1.247d+01,1.583d+03,3.935d+00,7.810d+00,6.558d-02,1.520d+00,1.084d-01)",
		"nej2" : "sigma_nej2 = sigma_v96(energy_eV,7.753d-01,5.708d+00,6.725d+01,1.005d+01,4.633d-01,7.654d+01,2.023d+00)",
		"nej3" : "sigma_nej3 = sigma_v96(energy_eV,5.566d+00,1.685d+03,6.409d+02,3.056d+00,8.290d-03,5.149d+00,6.687d+00)",
		"nej4" : "sigma_nej4 = sigma_v96(energy_eV,1.248d+00,2.430d+00,1.066d+02,8.999d+00,6.855d-01,9.169d+01,3.702d-01)",
		"nej5" : "sigma_nej5 = sigma_v96(energy_eV,1.499d+00,9.854d-01,1.350d+02,8.836d+00,1.656d+00,1.042d+02,1.435d+00)",
		"nej6" : "sigma_nej6 = sigma_v96(energy_eV,4.888d+00,1.198d+04,1.788d+00,2.550d+01,2.811d+01,2.536d-02,4.417d-02)",
		"nej7" : "sigma_nej7 = sigma_v96(energy_eV,1.003d+01,5.631d+01,3.628d+01,5.585d+00,0.000d+00,0.000d+00,0.000d+00)",
		"nej8" : "sigma_nej8 = sigma_v96(energy_eV,1.586d+02,6.695d+01,3.352d+01,2.002d+00,0.000d+00,0.000d+00,0.000d+00)",
		"nej9" : "sigma_nej9 = sigma_v96(energy_eV,4.304d+01,5.475d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"na" : "sigma_na = sigma_v96(energy_eV,6.139d+00,1.601d+00,6.148d+03,3.839d+00,0.000d+00,0.000d+00,0.000d+00)",
		"naj" : "sigma_naj = sigma_v96(energy_eV,8.203d+00,1.040d+03,8.259d+00,7.362d+00,2.328d+00,3.375d+00,4.010d+00)",
		"naj2" : "sigma_naj2 = sigma_v96(energy_eV,1.069d+01,1.885d+03,3.613d+00,9.803d+00,8.579d-02,3.725d+00,2.279d-01)",
		"naj3" : "sigma_naj3 = sigma_v96(energy_eV,6.690d-01,2.330d+00,1.205d+02,9.714d+00,7.365d-01,1.383d+02,4.260d+00)",
		"naj4" : "sigma_naj4 = sigma_v96(energy_eV,5.408d+00,2.346d+01,2.913d+01,8.260d+00,9.275d-01,2.204d+01,7.577d-01)",
		"naj5" : "sigma_naj5 = sigma_v96(energy_eV,4.846d+01,7.101d+01,3.945d+01,2.832d+00,1.285d-02,9.603d-04,6.378d-03)",
		"naj6" : "sigma_naj6 = sigma_v96(energy_eV,2.096d+00,1.609d+00,2.473d+02,7.681d+00,1.895d+00,9.940d+01,3.278d+00)",
		"naj7" : "sigma_naj7 = sigma_v96(energy_eV,1.535d+02,7.215d+03,3.886d-01,8.476d+00,9.121d-01,1.667d-01,1.766d-02)",
		"naj8" : "sigma_naj8 = sigma_v96(energy_eV,1.391d+01,4.729d+01,3.889d+01,5.265d+00,0.000d+00,0.000d+00,0.000d+00)",
		"naj9" : "sigma_naj9 = sigma_v96(energy_eV,2.268d+02,3.995d+01,5.315d+01,1.678d+00,0.000d+00,0.000d+00,0.000d+00)",
		"naj10" : "sigma_naj10 = sigma_v96(energy_eV,5.211d+01,4.525d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"mg" : "sigma_mg = sigma_v96(energy_eV,1.197d+01,1.372d+08,2.228d-01,1.574d+01,2.805d-01,0.000d+00,0.000d+00)",
		"mgj" : "sigma_mgj = sigma_v96(energy_eV,8.139d+00,3.278d+00,4.341d+07,3.610d+00,0.000d+00,0.000d+00,0.000d+00)",
		"mgj2" : "sigma_mgj2 = sigma_v96(energy_eV,1.086d+01,5.377d+02,9.779d+00,7.117d+00,2.604d+00,4.860d+00,3.722d+00)",
		"mgj3" : "sigma_mgj3 = sigma_v96(energy_eV,2.912d+01,1.394d+03,2.895d+00,6.487d+00,4.326d-02,9.402d-01,1.135d-01)",
		"mgj4" : "sigma_mgj4 = sigma_v96(energy_eV,9.762d-01,1.728d+00,9.184d+01,1.006d+01,8.090d-01,1.276d+02,3.979d+00)",
		"mgj5" : "sigma_mgj5 = sigma_v96(energy_eV,1.711d+00,2.185d+00,9.350d+01,9.202d+00,6.325d-01,1.007d+02,1.729d+00)",
		"mgj6" : "sigma_mgj6 = sigma_v96(energy_eV,3.570d+00,3.104d+00,6.060d+01,8.857d+00,1.422d+00,5.452d+01,2.078d+00)",
		"mgj7" : "sigma_mgj7 = sigma_v96(energy_eV,4.884d-01,6.344d-02,5.085d+02,9.385d+00,6.666d-01,5.348d+02,3.997d-03)",
		"mgj8" : "sigma_mgj8 = sigma_v96(energy_eV,3.482d+01,9.008d+02,1.823d+00,1.444d+01,2.751d+00,5.444d+00,7.918d-02)",
		"mgj9" : "sigma_mgj9 = sigma_v96(energy_eV,1.452d+01,4.427d+01,3.826d+01,5.460d+00,0.000d+00,0.000d+00,0.000d+00)",
		"mgj10" : "sigma_mgj10 = sigma_v96(energy_eV,2.042d+02,6.140d+01,2.778d+01,2.161d+00,0.000d+00,0.000d+00,0.000d+00)",
		"mgj11" : "sigma_mgj11 = sigma_v96(energy_eV,6.203d+01,3.802d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"al" : "sigma_al = sigma_v96(energy_eV,1.381d+01,7.195d+00,1.621d+03,3.642d+00,3.166d-01,2.041d-01,4.753d-01)",
		"alj" : "sigma_alj = sigma_v96(energy_eV,2.048d-01,6.948d-02,5.675d+02,9.049d+00,4.615d-01,9.149d+01,6.565d-01)",
		"alj2" : "sigma_alj2 = sigma_v96(energy_eV,1.027d+01,4.915d+00,1.990d+06,3.477d+00,0.000d+00,0.000d+00,0.000d+00)",
		"alj3" : "sigma_alj3 = sigma_v96(energy_eV,3.130d+00,1.513d+01,1.674d+01,1.180d+01,5.342d+00,3.994d+01,4.803d+00)",
		"alj4" : "sigma_alj4 = sigma_v96(energy_eV,2.414d+01,2.925d+02,6.973d+00,6.724d+00,1.000d-01,3.495d+00,2.701d-01)",
		"alj5" : "sigma_alj5 = sigma_v96(energy_eV,3.483d-01,1.962d-02,1.856d+01,2.084d+01,8.839d+00,5.675d-02,2.768d-01)",
		"alj6" : "sigma_alj6 = sigma_v96(energy_eV,2.636d+00,1.889d+02,1.338d+02,6.204d+00,1.836d+00,3.552d+01,8.223d-03)",
		"alj7" : "sigma_alj7 = sigma_v96(energy_eV,4.866d-01,2.350d-01,7.216d+02,8.659d+00,2.773d-01,5.704d+02,1.580d-01)",
		"alj8" : "sigma_alj8 = sigma_v96(energy_eV,1.842d+00,4.982d-01,2.568d+02,8.406d+00,6.945d-01,1.719d+02,6.595d+00)",
		"alj9" : "sigma_alj9 = sigma_v96(energy_eV,8.044d+00,1.774d+04,1.653d+00,2.655d+01,2.953d+01,2.538d-02,1.203d-02)",
		"alj10" : "sigma_alj10 = sigma_v96(energy_eV,2.355d+01,3.388d+01,3.432d+01,5.085d+00,0.000d+00,0.000d+00,0.000d+00)",
		"alj11" : "sigma_alj11 = sigma_v96(energy_eV,2.738d+02,4.036d+01,3.567d+01,1.915d+00,0.000d+00,0.000d+00,0.000d+00)",
		"alj12" : "sigma_alj12 = sigma_v96(energy_eV,7.281d+01,3.239d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"si" : "sigma_si = sigma_v96(energy_eV,2.317d+01,2.506d+01,2.057d+01,3.546d+00,2.837d-01,1.672d-05,4.207d-01)",
		"sij" : "sigma_sij = sigma_v96(energy_eV,2.556d+00,4.140d+00,1.337d+01,1.191d+01,1.570d+00,6.634d+00,1.272d-01)",
		"sij2" : "sigma_sij2 = sigma_v96(energy_eV,1.659d-01,5.790d-04,1.474d+02,1.336d+01,8.626d-01,9.613d+01,6.442d-01)",
		"sij3" : "sigma_sij3 = sigma_v96(energy_eV,1.288d+01,6.083d+00,1.356d+06,3.353d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sij4" : "sigma_sij4 = sigma_v96(energy_eV,7.761d-01,8.863d-01,1.541d+02,9.980d+00,1.303d+00,2.009d+02,4.537d+00)",
		"sij5" : "sigma_sij5 = sigma_v96(energy_eV,6.305d+01,7.293d+01,1.558d+02,2.400d+00,2.989d-03,1.115d+00,8.051d-02)",
		"sij6" : "sigma_sij6 = sigma_v96(energy_eV,3.277d-01,6.680d-02,4.132d+01,1.606d+01,3.280d+00,1.149d-02,6.396d-01)",
		"sij7" : "sigma_sij7 = sigma_v96(energy_eV,7.655d-01,3.477d-01,3.733d+02,8.986d+00,1.476d-03,3.850d+02,8.999d-02)",
		"sij8" : "sigma_sij8 = sigma_v96(energy_eV,3.343d-01,1.465d-01,1.404d+03,8.503d+00,1.646d+00,1.036d+03,2.936d-01)",
		"sij9" : "sigma_sij9 = sigma_v96(energy_eV,8.787d-01,1.950d-01,7.461d+02,8.302d+00,4.489d-01,4.528d+02,1.015d+00)",
		"sij10" : "sigma_sij10 = sigma_v96(energy_eV,1.205d+01,1.992d+04,1.582d+00,2.425d+01,2.392d+01,1.990d-02,1.007d-02)",
		"sij11" : "sigma_sij11 = sigma_v96(energy_eV,3.560d+01,2.539d+01,3.307d+01,4.728d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sij12" : "sigma_sij12 = sigma_v96(energy_eV,2.752d+02,4.754d+01,2.848d+01,2.135d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sij13" : "sigma_sij13 = sigma_v96(energy_eV,8.447d+01,2.793d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"s" : "sigma_s = sigma_v96(energy_eV,1.808d+01,4.564d+04,1.000d+00,1.361d+01,6.385d-01,9.935d-01,2.486d-01)",
		"sj" : "sigma_sj = sigma_v96(energy_eV,8.787d+00,3.136d+02,3.442d+00,1.281d+01,7.354d-01,2.782d+00,1.788d-01)",
		"sj2" : "sigma_sj2 = sigma_v96(energy_eV,2.027d+00,6.666d+00,5.454d+01,8.611d+00,4.109d+00,1.568d+01,9.421d+00)",
		"sj3" : "sigma_sj3 = sigma_v96(energy_eV,2.173d+00,2.606d+00,6.641d+01,8.655d+00,1.863d+00,1.975d+01,3.361d+00)",
		"sj4" : "sigma_sj4 = sigma_v96(energy_eV,1.713d-01,5.072d-04,1.986d+02,1.307d+01,7.880d-01,9.424d+01,6.265d-01)",
		"sj5" : "sigma_sj5 = sigma_v96(energy_eV,1.413d+01,9.139d+00,1.656d+03,3.626d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sj6" : "sigma_sj6 = sigma_v96(energy_eV,3.757d-01,5.703d-01,1.460d+02,1.135d+01,1.503d+00,2.222d+02,4.606d+00)",
		"sj7" : "sigma_sj7 = sigma_v96(energy_eV,1.462d+01,3.161d+01,1.611d+01,8.642d+00,1.153d-03,1.869d+01,3.037d-01)",
		"sj8" : "sigma_sj8 = sigma_v96(energy_eV,1.526d-01,9.646d+03,1.438d+03,5.977d+00,1.492d+00,1.615d-03,4.049d-01)",
		"sj9" : "sigma_sj9 = sigma_v96(energy_eV,1.040d+01,5.364d+01,3.641d+01,7.090d+00,2.310d+00,1.775d+01,1.663d+00)",
		"sj10" : "sigma_sj10 = sigma_v96(energy_eV,6.485d+00,1.275d+01,6.583d+01,7.692d+00,1.678d+00,3.426d+01,1.370d-01)",
		"sj11" : "sigma_sj11 = sigma_v96(energy_eV,2.443d+00,3.490d-01,5.411d+02,7.769d+00,7.033d-01,2.279d+02,1.172d+00)",
		"sj12" : "sigma_sj12 = sigma_v96(energy_eV,1.474d+01,2.294d+04,1.529d+00,2.568d+01,2.738d+01,2.203d-02,1.073d-02)",
		"sj13" : "sigma_sj13 = sigma_v96(energy_eV,3.310d+01,2.555d+01,3.821d+01,5.037d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sj14" : "sigma_sj14 = sigma_v96(energy_eV,4.390d+02,2.453d+01,4.405d+01,1.765d+00,0.000d+00,0.000d+00,0.000d+00)",
		"sj15" : "sigma_sj15 = sigma_v96(energy_eV,1.104d+02,2.139d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"ar" : "sigma_ar = sigma_v96(energy_eV,1.709d+01,2.106d+01,2.645d+02,4.796d+00,4.185d-01,1.688d+00,8.943d-01)",
		"arj" : "sigma_arj = sigma_v96(energy_eV,2.494d+01,2.503d+01,1.272d+02,4.288d+00,5.108d-01,9.299d-01,7.195d-01)",
		"arj2" : "sigma_arj2 = sigma_v96(energy_eV,1.417d+01,3.580d+01,3.776d+01,5.742d+00,6.316d-01,2.384d+00,1.794d+00)",
		"arj3" : "sigma_arj3 = sigma_v96(energy_eV,6.953d+00,2.035d+01,1.400d+01,9.595d+00,8.842d-01,7.501d+00,1.806d-01)",
		"arj4" : "sigma_arj4 = sigma_v96(energy_eV,1.031d+01,9.946d+00,7.444d+01,6.261d+00,4.885d-01,6.406d+00,3.659d-03)",
		"arj5" : "sigma_arj5 = sigma_v96(energy_eV,5.440d-01,1.080d+00,9.419d+02,7.582d+00,1.107d+01,1.700d+02,1.587d+01)",
		"arj6" : "sigma_arj6 = sigma_v96(energy_eV,2.966d-02,3.693d+00,9.951d+03,7.313d+00,1.363d-02,4.383d-04,2.513d+00)",
		"arj7" : "sigma_arj7 = sigma_v96(energy_eV,3.884d+00,3.295d+01,7.082d+02,4.645d+00,0.000d+00,0.000d+00,0.000d+00)",
		"arj8" : "sigma_arj8 = sigma_v96(energy_eV,1.926d-01,8.279d-01,2.392d+02,1.121d+01,1.434d+00,3.814d+01,4.649d+00)",
		"arj9" : "sigma_arj9 = sigma_v96(energy_eV,1.040d+01,8.204d+00,1.495d+01,1.115d+01,9.203d-04,3.804d+01,6.390d-01)",
		"arj10" : "sigma_arj10 = sigma_v96(energy_eV,1.257d-01,1.760d+03,1.579d+03,6.714d+00,1.975d+00,3.286d-03,3.226d-01)",
		"arj11" : "sigma_arj11 = sigma_v96(energy_eV,5.310d+00,7.018d-01,1.001d+02,8.939d+00,4.987d-01,1.099d+02,2.202d-01)",
		"arj12" : "sigma_arj12 = sigma_v96(energy_eV,3.209d-01,2.459d-02,2.285d+03,8.810d+00,6.692d-01,2.068d+03,2.113d+01)",
		"arj13" : "sigma_arj13 = sigma_v96(energy_eV,1.557d+00,4.997d-02,5.031d+02,8.966d+00,2.938d-01,4.552d+02,6.459d+00)",
		"arj14" : "sigma_arj14 = sigma_v96(energy_eV,1.888d+01,2.571d+04,1.475d+00,2.634d+01,2.909d+01,2.445d-02,1.054d-02)",
		"arj15" : "sigma_arj15 = sigma_v96(energy_eV,4.154d+01,2.135d+01,4.118d+01,4.945d+00,0.000d+00,0.000d+00,0.000d+00)",
		"arj16" : "sigma_arj16 = sigma_v96(energy_eV,4.468d+02,3.108d+01,3.039d+01,2.092d+00,0.000d+00,0.000d+00,0.000d+00)",
		"arj17" : "sigma_arj17 = sigma_v96(energy_eV,1.399d+02,1.690d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"ca" : "sigma_ca = sigma_v96(energy_eV,1.278d+01,5.370d+05,3.162d-01,1.242d+01,4.477d-01,1.012d-03,1.851d-02)",
		"caj" : "sigma_caj = sigma_v96(energy_eV,1.553d+01,1.064d+07,7.790d-01,2.130d+01,6.453d-01,2.161d-03,6.706d-02)",
		"caj2" : "sigma_caj2 = sigma_v96(energy_eV,2.436d+01,3.815d+01,2.931d+02,3.944d+00,3.126d-01,1.802d+00,1.233d+00)",
		"caj3" : "sigma_caj3 = sigma_v96(energy_eV,4.255d+00,7.736d+00,1.355d+01,1.236d+01,1.369d+00,1.467d+01,3.298d-02)",
		"caj4" : "sigma_caj4 = sigma_v96(energy_eV,6.882d-01,1.523d-01,1.502d+02,1.061d+01,8.227d+00,1.210d+02,3.876d+00)",
		"caj5" : "sigma_caj5 = sigma_v96(energy_eV,9.515d+00,7.642d+01,8.973d+01,5.141d+00,2.471d+00,4.829d+00,5.824d+00)",
		"caj6" : "sigma_caj6 = sigma_v96(energy_eV,8.080d-01,4.760d-01,3.682d+02,8.634d+00,5.720d-01,1.487d+02,1.283d+00)",
		"caj7" : "sigma_caj7 = sigma_v96(energy_eV,1.366d+00,6.641d-01,3.188d+02,8.138d+00,2.806d-01,1.039d+02,3.329d+00)",
		"caj8" : "sigma_caj8 = sigma_v96(energy_eV,5.520d-02,2.076d+02,1.790d+04,5.893d+00,1.843d-03,2.826d-04,1.657d+00)",
		"caj9" : "sigma_caj9 = sigma_v96(energy_eV,1.605d+01,1.437d+01,6.989d+02,3.857d+00,0.000d+00,0.000d+00,0.000d+00)",
		"caj10" : "sigma_caj10 = sigma_v96(energy_eV,2.288d-01,9.384d-01,2.549d+02,1.103d+01,1.390d+00,2.478d+01,3.100d+00)",
		"caj11" : "sigma_caj11 = sigma_v96(energy_eV,2.345d+01,1.227d+01,1.312d+01,9.771d+00,6.842d-04,2.417d+01,5.469d-01)",
		"caj12" : "sigma_caj12 = sigma_v96(energy_eV,1.008d+01,1.849d+03,1.792d+04,2.868d+00,2.410d+02,6.138d-03,6.931d+01)",
		"caj13" : "sigma_caj13 = sigma_v96(energy_eV,9.980d+00,1.116d+00,5.918d+01,9.005d+00,3.879d+00,7.104d+01,5.311d+00)",
		"caj14" : "sigma_caj14 = sigma_v96(energy_eV,1.309d+02,5.513d+01,3.828d+02,2.023d+00,9.084d-02,1.833d-02,9.359d-01)",
		"caj15" : "sigma_caj15 = sigma_v96(energy_eV,4.293d+00,1.293d+00,1.691d+01,1.438d+01,3.461d-05,9.363d-01,4.589d-02)",
		"caj16" : "sigma_caj16 = sigma_v96(energy_eV,2.618d+01,2.028d+04,1.456d+00,2.560d+01,2.803d+01,2.402d-02,9.323d-03)",
		"caj17" : "sigma_caj17 = sigma_v96(energy_eV,9.472d+01,1.105d+01,3.818d+01,4.192d+00,0.000d+00,0.000d+00,0.000d+00)",
		"caj18" : "sigma_caj18 = sigma_v96(energy_eV,6.297d+02,1.936d+01,3.921d+01,1.862d+00,0.000d+00,0.000d+00,0.000d+00)",
		"caj19" : "sigma_caj19 = sigma_v96(energy_eV,1.729d+02,1.369d+02,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fe" : "sigma_fe = sigma_v96(energy_eV,5.461d-02,3.062d-01,2.671d+07,7.923d+00,2.069d+01,1.382d+02,2.481d-01)",
		"fej" : "sigma_fej = sigma_v96(energy_eV,1.761d-01,4.365d+03,6.298d+03,5.204d+00,1.141d+01,9.272d+01,1.075d+02)",
		"fej2" : "sigma_fej2 = sigma_v96(energy_eV,1.698d-01,6.107d+00,1.555d+03,8.055d+00,8.698d+00,1.760d+02,1.847d+01)",
		"fej3" : "sigma_fej3 = sigma_v96(energy_eV,2.544d+01,3.653d+02,8.913d+00,6.538d+00,5.602d-01,0.000d+00,0.000d+00)",
		"fej4" : "sigma_fej4 = sigma_v96(energy_eV,7.256d-01,1.523d-03,3.736d+01,1.767d+01,5.064d+01,8.871d+01,5.280d-02)",
		"fej5" : "sigma_fej5 = sigma_v96(energy_eV,2.656d+00,5.259d-01,1.450d+01,1.632d+01,1.558d+01,3.361d+01,3.743d-03)",
		"fej6" : "sigma_fej6 = sigma_v96(energy_eV,5.059d+00,2.420d+04,4.850d+04,2.374d+00,2.516d-03,4.546d-01,2.683d+01)",
		"fej7" : "sigma_fej7 = sigma_v96(energy_eV,7.098d-02,1.979d+01,1.745d+04,6.750d+00,2.158d+02,2.542d+03,4.672d+02)",
		"fej8" : "sigma_fej8 = sigma_v96(energy_eV,6.741d+00,2.687d+01,1.807d+02,6.290d+00,2.387d-04,2.494d+01,8.251d+00)",
		"fej9" : "sigma_fej9 = sigma_v96(energy_eV,6.886d+01,6.470d+01,2.062d+01,4.111d+00,2.778d-04,1.190d-05,6.570d-03)",
		"fej10" : "sigma_fej10 = sigma_v96(energy_eV,8.284d+00,3.281d+00,5.360d+01,8.571d+00,3.279d-01,2.971d+01,5.220d-01)",
		"fej11" : "sigma_fej11 = sigma_v96(energy_eV,6.295d+00,1.738d+00,1.130d+02,8.037d+00,3.096d-01,4.671d+01,1.425d-01)",
		"fej12" : "sigma_fej12 = sigma_v96(energy_eV,1.317d-01,2.791d-03,2.487d+03,9.791d+00,6.938d-01,2.170d+03,6.852d-03)",
		"fej13" : "sigma_fej13 = sigma_v96(energy_eV,8.509d-01,1.454d-01,1.239d+03,8.066d+00,4.937d-01,4.505d+02,2.504d+00)",
		"fej14" : "sigma_fej14 = sigma_v96(energy_eV,5.555d-02,2.108d+02,2.045d+04,6.033d+00,1.885d-03,2.706d-04,1.628d+00)",
		"fej15" : "sigma_fej15 = sigma_v96(energy_eV,2.873d+01,1.207d+01,5.150d+02,3.846d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fej16" : "sigma_fej16 = sigma_v96(energy_eV,3.444d-01,1.452d+00,3.960d+02,1.013d+01,1.264d+00,2.891d+01,3.404d+00)",
		"fej17" : "sigma_fej17 = sigma_v96(energy_eV,3.190d+01,2.388d+00,2.186d+01,9.589d+00,2.902d-02,3.805d+01,4.805d-01)",
		"fej18" : "sigma_fej18 = sigma_v96(energy_eV,7.519d-04,6.066d-05,1.606d+06,8.813d+00,4.398d+00,1.915d+06,3.140d+01)",
		"fej19" : "sigma_fej19 = sigma_v96(energy_eV,2.011d+01,4.455d-01,4.236d+01,9.724d+00,2.757d+00,6.847d+01,3.989d+00)",
		"fej20" : "sigma_fej20 = sigma_v96(energy_eV,9.243d+00,1.098d+01,7.637d+01,7.962d+00,1.748d+00,4.446d+01,3.512d+00)",
		"fej21" : "sigma_fej21 = sigma_v96(energy_eV,9.713d+00,7.204d-02,1.853d+02,8.843d+00,9.551d-03,1.702d+02,4.263d+00)",
		"fej22" : "sigma_fej22 = sigma_v96(energy_eV,4.575d+01,2.580d+04,1.358d+00,2.604d+01,2.723d+01,3.582d-02,8.712d-03)",
		"fej23" : "sigma_fej23 = sigma_v96(energy_eV,7.326d+01,1.276d+01,4.914d+01,4.941d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fej24" : "sigma_fej24 = sigma_v96(energy_eV,1.057d+03,1.195d+01,5.769d+01,1.718d+00,0.000d+00,0.000d+00,0.000d+00)",
		"fej25" : "sigma_fej25 = sigma_v96(energy_eV,2.932d+02,8.099d+01,3.288d+01,2.963d+00,0.000d+00,0.000d+00,0.000d+00)"
	}
	if(myatom in cross): return cross[myatom]
	return "sigma_" + myatom + " = 0.d0"

	
################################
#create photochemistry parts
# for a given reaction (react)
def get_ph_stuff(react):
	if(len(react.reactants)>1): return None #die("ERROR:too much reactants for a photoreaction!", react.verbatim)
	if(react.kphrate==None): return None
	myreag = react.reactants
	mol = myreag[0]
	if(mol.charge<0): return None #negative species do not partecipate to photoheating
	if(mol.charge>=0 and mol.is_atom):
		reaname = mol.phname.capitalize()
		qromos = get_photo_qromosV96(mol.phname)
		heat = get_photo_heatV96(mol.phname)
		cross = get_photo_crossV96(mol.phname)
	else:
		reaname = mol.phname.capitalize() + "R" + str(react.idx)
		Emin = react.Tmin
		Emax = react.Tmax
		qromos = "call qromos(intf, sigma_"+reaname+", "+str(Emin)+", "+str(Emax)+", krome_kph_"+reaname+", midsqls)"
		heat = "heat_"+reaname+" = ("+react.kphrate+") * (energy_eV - "+str(Emin)+")"
		cross = "sigma_"+reaname+" = "+react.kphrate
	#creates cross section function, e.g. sigma_H(energy_eV)
	ph_func = ""
	ph_func += "!************************\n"
	ph_func += "function sigma_"+reaname+"(energy_eV)\n"
	ph_func += "\treal*8::sigma_"+reaname+",energy_eV\n"
	ph_func += "\t" + cross + "\n"
	ph_func += "end function sigma_"+reaname+"\n"
	ph_func += "\n"
	
	#creates photoheating function, e.g. heat_H(energy_eV)
	ph_heat = ""
	ph_heat += "!************************\n"
	ph_heat += "function heat_"+reaname+"(energy_eV)\n"
	ph_heat += "\treal*8::heat_"+reaname+",energy_eV\n"
	ph_heat += "\t" + heat + "\n"
	ph_heat += "end function heat_"+reaname+"\n"
	ph_heat += "\n"

	return {"reaname":reaname, "ph_func":ph_func, "ph_heat":ph_heat, "qromos":qromos}


#################################
# returns the licence of KROME
def get_licence_header(version, codename, short=False):
	import datetime
	header =  """!!*************************************************************
	!! This file has been generated with:
	!! krome #version# "#codename#" on #date#.
	!!
	!!KROME is a nice and friendly chemistry package for a wide range of 
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
	!!Written and developed by Tommaso Grassi
	!! tommasograssi@gmail.com,
	!! Starplan Center, Copenhagen.
	!! Niels Bohr Institute, Copenhagen.
	!!
	!!Co-developer Stefano Bovino
  	!! sbovino@astro.physik.uni-goettingen.de
	!! Institut fuer Astrophysik, Goettingen.
	!!
	!!Others (alphabetically): F.A. Gianturco, T. Haugboelle, 
	!! J.Prieto, D.R.G. Schleicher, D. Seifried, E. Simoncini, 
	!! E. Tognelli
	!!
	!!
	!!KROME is provided \"as it is\", without any warranty. 
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

	if(short): header = """!!*************************************************************
	!!This file has been generated with:
	!!KROME #version# on #date#
	!!see http://kromepackage.org
	!!
	!!Written and developed by Tommaso Grassi
	!!
	!!Co-developer Stefano Bovino
	!!Others (alphabetically): F.A. Gianturco, T. Haugboelle, 
	!! J.Prieto, D.R.G. Schleicher, D. Seifried, E. Simoncini, 
	!! E. Tognelli.
	!!KROME is provided \"as it is\", without any warranty.
	!!*************************************************************\n"""

	datenow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	header = header.replace("#date#",datenow).replace("#version#",version).replace("#codename#",codename)
	return header.replace("\t","").replace("!!","   ! ")

#################################
#breaks a string (mystr), in piece
# of length (sublen), using a 
# separator (sep)
def trunc(mystr,sublen,sep):
	astr = mystr.split(sep)
	s = z = ""
	for x in astr:
		z += x + sep
		s += x + sep
		if(len(z)>sublen):
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
		if(line[:len(arg)]==arg): return True
	return False

###########################
#this function returns if the string line
# ends with one of the items in the array 
# of string aarg
def lend(aarg,line):
	for arg in aarg:
		if(line[len(line)-len(arg):]==arg): return True
	return False

###############################
#this function indent f90 file and remove multiple blank lines
def indentF90(filename):
	import os
	#check if the file exists else return
	if(not(os.path.isfile(filename))): return
	
	#open file for indent
	fh = open(filename,"rb")
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
		if(module_head in srow): module_head_found = True #do not duplicate module header
		#check module begin
		if(lbeg(["module"], srow) and not(module_head_found)): 
			arow.append("\n") #blank line
			arow.append(module_head+"\n") #comment
		if(lbeg(tokenclose, srow)): nind -= 1 #check if the line ends with one of tokenclose
		indent = (" "*(nind*nspace)) #compute number of spaces for indent
		if(is_amper): indent = (" "*(2*nspace)) + indent #increas indent in case of previous &
		if(srow!=""):
			if(srow[0]=="#"): indent = "" #no indent for pragmas
		if(not(srow=="" and is_blank)): arow.append(indent+srow+"\n") #append indented line to array of rows
		is_amper = False #is a line after ampersend flag
		if(lend(["&"], srow)): is_amper = True #check if the line ends with &
		is_blank = (srow=="") #flag for blank line mode
		if(lbeg(tokenopen, srow)): nind += 1 #check if the line ends with one of tokenclose
		if(lbeg(["if"],srow) and "then" in srow): nind += 1 #check if line stats with if and has then
	fh.close()

	#write the new file
	fh = open(filename,"w")
	for x in arow:
		fh.write(x)
	fh.close()

################################
#This function writes an error 
# and exit
def die(msg):
	import sys
	print msg
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
	 + " create bigger and better idiots.  So far the Universe is winning.","Rich Cook"],
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
	["Chemistry has been termed by the physicist as the messy part of physics", "Frederick Soddy "]
	]
	qrange = 1
	print 
	if(qall): qrange = len(quotes)
	for i in range(qrange):
		irand = int(random.random()*(len(quotes)))
		if(qall): irand = i
		qtup = quotes[irand]
		myqt = trunc(str(irand+1)+". "+qtup[0],40," ").upper().strip()
		amyqt = myqt.split("\n")
		lqt = max([len(x) for x in amyqt])
		print
		if(i==0): print "*"*lqt
		print myqt
		if(qtup[1].strip()==""): qtup[1] = "Anonymous"
		print "--- "+qtup[1]
		if(i==qrange-1): 
			print "*"*lqt
			print

