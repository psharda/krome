import os,glob,shutil,argparse,re
from kromelib import *
from os import listdir
from os.path import isfile, join


class krome():
	#set defaults
	solver_MF = 222
	force_rwork = useHeating = doReport = checkConserv = useFileIdx = buildCompact = useEquilibrium = False
	use_implicit_RHS = use_photons = useTabs = useDvodeF90 = useTopology = useFlux = skipDup = False
	useCoolingAtomic = useCoolingH2 = useCoolingH2GP98 = useCoolingHD = useCoolingZ = use_cooling = useCoolingDust = useCoolingCont = False
	useCoolingCompton = useH2opacity = useCoolingCIE = False
	useCoolingZC = useCoolingZCp = useCoolingZSi = useCoolingZSip = useCoolingZO = useCoolingZOp = useCoolingZFe = useCoolingZFep = False
	useReverse = useCustomCoe = useODEConstant = cleanBuild = usePlainIsotopes = useDust = use_thermo = False
	usePhIoniz = useHeatingCompress = useHeatingPhoto = useHeatingChem = useDecoupled = useCoolingdH = useHeatingdH = useCoolingChem = False
	pedanticMakefile = useFakeOpacity = False
	useX = has_plot = doIndent = useTlimits = useODEthermo = True
	useDustGrowth = useDustSputter = useDustH2 = useDustT = False
	doRamses = doFlash = doEnzo = wrapC = False
	typeGamma = "DEFAULT"
	test_name = "default"
	is_test = False
	TlimitOpLow = "GE"
	TlimitOpHigh = "LT"
	customCoeFunction = "[CUSTOM COE FUNCTION NOT SET!]"
	buildFolder = "build/"
	TminAuto = 1e99
	TmaxAuto = 0e0
	RTOL = 1e-4 #default relative tolerance
	ATOL = 1e-20 #default absolute tolerance
	dustArraySize = dustTypesSize = 0
	maxord = 0
	dustTypes = []
	specs = []
	reacts = []
	constantList = []
	dummy = molec()
	coevars = dict() #variables in function coe() (krome_subs.f90)
	coevarsODE = dict() #variables in function fex() (krome_ode.f90)
	implicit_arrays = totMetals = ""
	thermodata = dict() #thermochemistry data (nasa polynomials)
	parser = filename = ""
	deltajacMode = "RELATIVE" #increment mode: RELATIVE or ABSOLUTE
	deltajac = "1d-3" #increment (relative or absolute, see deltajacMode)
	atols = [] #custom ATOLs
	rtols = [] #custom RTOLs
	customODEs = [] #custom ODEs
	version = "13.11"
	codename = "Astonishing Ansatz"
	#########################################
	def init_argparser(self):

		tests = ", ".join(os.walk('tests').next()[1])
		self.parser = argparse.ArgumentParser(description="KROME a package for astrochemistry")
		self.parser.add_argument("-n", help="reaction network file", metavar='FILENAME')
		self.parser.add_argument("-forceMF21", action="store_true", help="force explicit sparsity and Jacobian")
		self.parser.add_argument("-forceMF222", action="store_true", help="force internal-generated sparsity and Jacobian")
	 
		self.parser.add_argument("-test",help=("Create a test model in /build. TEST can be: "+tests+"."))
		self.parser.add_argument("-heating", metavar='TERMS', help="heating options, TERMS can be COMPRESS, PHOTO, CHEM, DH")
		self.parser.add_argument("-cooling", metavar='TERMS', help="cooling options, TERMS can be ATOMIC, H2, HD, Z, DH, DUST, H2GP98, COMPTON, CIE, CI, CII, SiI, SiII, OI, OII, FeI, FeII (e.g. -cooling ATOMIC,CII,OI,FeI), CHEM")
		self.parser.add_argument("-useN", action="store_true",help="use number densities (1/cm3) as input/ouput instead of fractions (#)")
		self.parser.add_argument("-useH2opacity", action="store_true",help="use H2 opacity for H2 cooling")
		self.parser.add_argument("-gamma",help="define the adiabatic index according to OPTION that can be FULL for employing Grassi et al. 2011, or a custom F90 expression e.g. -gamma 5.d0/3.d0",metavar="OPTION")
		self.parser.add_argument("-iRHS", action="store_true", help="implicit loop-based RHS (suggested for large systems)")

		self.parser.add_argument("-maxord", help="max order of the BDF solver. Default (and maximum values) is 5.")
		self.parser.add_argument("-customATOL", help="file with the list of the individual ATOLs in the form SPECIES ATOL in each line, e.g. H2 1d-20", metavar="filename")
		self.parser.add_argument("-customRTOL", help="file with the list of the individual RTOLs in the form SPECIES RTOL in each line, e.g. H3+ 1d-4", metavar="filename")
		self.parser.add_argument("-ATOL", help="set solver absolute tolerance to the float or double value ATOL, e.g. -atol 1d-40 Default is ATOL=1d-20")
		self.parser.add_argument("-RTOL", help="set solver relative tolerance to the float double value RTOL, e.g. -rtol 1e-5 Default is RTOL=1d-4")
		self.parser.add_argument("-usePhIoniz", action="store_true", help="include photochemistry")
		self.parser.add_argument("-useEquilibrium", action="store_true", help="check if the solver has reached the equilbirum. If so break the solver's loop and return the values found. It is useful when the system oscillates around a solution (as in some photoheating cases). To be used with caution.")
		self.parser.add_argument("-dust", help="include dust ODE using N bins for each TYPE, e.g. -dust 10,C,Si set 10 dust carbon bins and 10 dust silicon dust bins. Require a call to the krome_init_dust subroutine. See test=dust for an example.")
		self.parser.add_argument("-dustOptions", help="activate dust options: (GROWTH) dust growth, (SPUTTER) sputtering, (H2) molecular hydrogen formation on dust, and (T) dust temperature. The last option provide a template for the FEX routine.", metavar="OPTIONS")
		self.parser.add_argument("-compact", action="store_true", help="creates a single fortran file with all the modules instead of various file with the different modules. Solver files remain stand-alone (see example make in test/MakefileCompact).")
		self.parser.add_argument("-useDvodeF90", action="store_true", help="use Dvode implementation in F90 (slower)")
		self.parser.add_argument("-useTabs", action="store_true", help="use tabulated rate coefficients (free parameter: temperature)")
		self.parser.add_argument("-report", action="store_true", help="generate report file in the main call to krome as KROME_ERROR_REPORT and when calling the fex as KROME_ODE_REPORT. It also stores abundances evolution in fex as fort.98, and prepares a report.gps gnuplot script file to plot evolutions callable in gnuplot with load 'report.gps'. Warning: it slows the whole system!")
		self.parser.add_argument("-checkConserv", action="store_true", help="check mass conservation during integration (slower)")
		self.parser.add_argument("-useFileIdx", action="store_true", help="use the reaction index in the reaction file instead of using the automatic progressive index starting from 1. Useful with rate coefficients that depends on other coefficients, e.g. k(10) = 1d-2*k(3)")
		self.parser.add_argument("-skipDup", action="store_true", help="skip duplicate reactions")
		self.parser.add_argument("-Tlimit", metavar="opLow,opHigh", help="set the operators for all the reaction temperature limits where opLow is the operator for the first temperature value in the reaction file, and opHigh is for the second one. e.g. if the T limits for a 	given reaction are 10. and 1d4 the option -Tlmit GE,LE will provide (Tgas>=10. AND Tgas<=1d4) as the reaction range of validity. Operators opLow and opHigh must be one of the following: LE, GE, LT, GT.")
		self.parser.add_argument("-noTlimits", action="store_true", help="ignore rate coefficient temperature limits.")
		self.parser.add_argument("-reverse", action="store_true", help="create reverse reaction from the current network using NASA polynomials.")
		self.parser.add_argument("-useCustomCoe", help="use a user-defined custom function that returns a real*8 array of size NREA = number of reactions, that replaces the standard rate coefficient calculation function. Note that FUNCTION must be explicitly included in krome_user_commons module.", metavar="FUNCTION")
		self.parser.add_argument("-useODEConstant", help="postpone an expression to each ODE. EXPRESSION must be a valid f90 expression (e.g. *3.d0 or +1.d-10)", metavar="EXPRESSION")
		self.parser.add_argument("-skipODEthermo", action="store_true", help="do not compute dT/dt in the ODE RHS function (fex)")
		self.parser.add_argument("-usePlainIsotopes", action="store_true", help="use kA format for isotopes instead of [k]A format, where k is the isotopic number and A is the atom name, e.g. krome looks for 14C instead of [14]C in the reactions file.")
		self.parser.add_argument("-project", help="build everything in a folder called build_NAME instead of building all in the default build folder. It also creates a NAME.kpj file with the krome input used.",metavar="NAME")
		self.parser.add_argument("-clean", action="store_true", help="clean all in /build (including krome_user_commons.f90 that is normally kept by default) before creating new f90 files.")
		self.parser.add_argument("-pedantic", action="store_true", help="uses a pedantic Makefile (debug purposes)")
		self.parser.add_argument("-forceRWORK", help="force the size of RWORK to N", metavar="N")
		self.parser.add_argument("-ramses", action="store_true", help="create patches for RAMSES")
		self.parser.add_argument("-flash", action="store_true", help="create patches for FLASH")
		self.parser.add_argument("-enzo", action="store_true", help="create patches for ENZO")
		self.parser.add_argument("-C", action="store_true", help="create a simple C wrapper")
		self.parser.add_argument("-customODE", help="file with the list of custom ODEs", metavar="FILENAME")
		
	
	######################################
	#select test name
	def select_test(self,argv):
		parser = self.parser
		args = parser.parse_args()

		if(args.test):
			self.is_test = True
		else:
			return
		#test_name = (arg.strip().replace("-test=",""))
		#print "Reading option -test (test="+test_name+")"
		if(args.test=="cloud"):
			[argv.append(x) for x in ["-useN","-iRHS"]]
			filename = "networks/react_cloud"
		elif(args.test=="slowmanifold"):
			[argv.append(x) for x in ["-useN"]]
			filename = "networks/react_SM"
		elif(args.test=="shock1Dcool"):
			[argv.append(x) for x in ["-cooling=H2,HD,Z,DH"]]
			filename = "networks/react_primordial"
		elif(args.test=="shock1D"):
			filename = "networks/react_primordial"
		elif(args.test=="shock1Dphoto"):
			[argv.append(x) for x in ["-usePhIoniz","-heating=PHOTO","-cooling=ATOMIC,H2,HD,Z","-useEquilibrium"]]
			filename = "networks/react_primordial_photo"
		elif(args.test=="shock1Dlarge"):
			[argv.append(x) for x in ["-iRHS"]]
			filename = "networks/react_WH2008"
		elif(args.test=="dust"):
			[argv.append(x) for x in ["-dust=10,C,Si","-useN","-dustOptions=GROWTH,SPUTTER"]]
			filename = "networks/react_primordial"
		elif(args.test=="compact"):
			[argv.append(x) for x in ["-compact"]]
			filename = "networks/react_primordial"
		elif(args.test=="map"):
			[argv.append(x) for x in ["-cooling=ATOMIC,HD,H2", "-heating=PHOTO","-usePhIoniz"]]
			filename = "networks/react_primordial_photoH2"
		elif(args.test=="collapse"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CONT,CHEM", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL"]]
			filename = "networks/react_primordial2"
		elif(args.test=="collapseZ"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CI,CII,OI,OII,SiII,FeII,CONT,CHEM", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL","-ATOL=1d-40","-maxord=1"]]
			filename = "networks/react_primordialZ2"
		elif(args.test=="collapseUV"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CIE,ATOMIC", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-useN","-gamma=FULL"]]
			filename = "networks/react_primordial_UV"
		elif(args.test=="collapseDUST"):
			[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE,DUST,HD", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL","-dust=1,C","-dustOptions=H2"]]
			filename = "networks/react_primordial"
		elif(args.test=="reverse"):
			[argv.append(x) for x in ["-useN","-reverse"]]
			filename = "networks/react_NO"
		elif(args.test=="atmosphere"):
			[argv.append(x) for x in ["-useN"]]
			filename = "networks/react_kast80"
		elif(args.test=="wrapC"):
			[argv.append(x) for x in ["-useN"]]
			filename = "networks/react_primordial2"
		elif(args.test=="lotkav"):
			[argv.append(x) for x in ["-useN","-customODE=tests/lotkav/lotkav"]]
			filename = "networks/react_dummy"
		else:
			tests = ", ".join(os.walk('tests').next()[1])
			print "ERROR: test \""+args.test+"\" not present!"
			print "Available tests are: "+tests
			sys.exit()
		self.filename = filename
		self.test_name = args.test

	##########################################
	def argparsing(self,argv):
		args = self.parser.parse_args()


		#you can select only one -forceMF
		if((args.forceMF222) and (args.forceMF21)):
			die("ERROR: options -forceMF222 and -forceMF21 are mutually exclusive: choose one.")

		#get filename
		if(not(self.is_test) and args.n): self.filename = args.n
		#chech if reactions file exists
		if(not(os.path.isfile(self.filename))): die("ERROR: Reaction file \""+self.filename+"\" doesn't exist!")


		#set implicit RHS
		if(args.iRHS):
			self.use_implicit_RHS = True
			self.solver_MF = 222
			print "Reading option -iRHS"
		#force MF=21
		if(args.forceMF21):
			self.solver_MF = 21
			print "Reading option -forceMF21"
		#force MF=222
		if(args.forceMF222):
			self.solver_MF = 222
			print "Reading option -forceMF222"
		#use numeric density instead of fractions as input
		if(args.useN):
			self.useX = False
			print "Reading option -useN"
		#force to use photons
		#if(args.usePhot):
		#	self.use_photons = True
		#	print "Reading option -usePhot"
		#use rate tables
		if(args.useTabs):
			self.useTabs = True
			print "Reading option -useTabs"
		#use f90 solver
		if(args.useDvodeF90):
			self.useDvodeF90 = True
			print "Reading option -useDvodeF90"
		#do report
		if(args.report):
			self.doReport = True
			print "Reading option -report"
		#check mass conservation
		if(args.checkConserv):
			self.checkConserv = True
			print "Reading option -checkConserv"
		#use reaction indexes in reaction file
		if(args.useFileIdx):
			self.useFileIdx = True
			print "Reading option -useFileIdx"
		#write a single compact file krome_all.f90
		if(args.compact):
			self.buildCompact = True
			print "Reading option -compact"
		#perform a clean build
		if(args.clean):
			self.cleanBuild = True
			print "Reading option -clean"
		#build isotopes automatically
		if(args.usePlainIsotopes):
			self.usePlainIsotopes = True
			print "Reading option -usePlainIsotopes"
		#use photoionization from Verner et al. 1996
		if(args.usePhIoniz):
			self.usePhIoniz = True
			print "Reading option -usePhIoniz"
		#use equilibrium check to break loops earlier
		if(args.useEquilibrium):
			self.useEquilibrium = True
			print "Reading option -useEquilibrium"
		#do not use temperature limits
		if(args.noTlimits):
			self.useTlimits = False
			print "Reading option -noTlimits"
		#skip duplicated reactions
		if(args.skipDup):
			self.skipDup = True
			print "Reading option -skipDup"
		#skip duplicated reactions
		if(args.pedantic):
			self.pedanticMakefile = True
			print "Reading option -pedantic"
		#use reverse kinetics
		if(args.reverse):
			self.useReverse = True
			print "Reading option -reverse"
		#use H2opacity following 
		if(args.useH2opacity):
			self.useH2opacity = True
			print "Reading option -useH2opacity"

		#use cooling dT/dt in the ODE fex
		if(args.skipODEthermo):
			self.useODEthermo = False
			print "Reading option -skipODEthermo"

		#creates ramses patches
		if(args.ramses):
			self.doRamses = True
			print "Reading option -ramses"
			if(not(args.compact)):
				print "ERROR: the patch for RAMSES requires the -compact option!"
				sys.exit()
			if(self.is_test):
				print "ERROR: -test option and -ramses are incompatible!"
				sys.exit()
			if(self.useX):
				print "ERROR: the patch for RAMSES requires the -useN option!"
				sys.exit()
		#creates flash patches
		if(args.flash):
			self.doFlash = True
			print "Reading option -flash"

		#creates enzo patches
		if(args.enzo):
			self.doEnzo = True
			print "Reading option -enzo"

		#creates simple C wrapper
		if(args.C):
			self.wrapC = True
			print "Reading option -C"


		#determine Tgas limit operators
		if(args.Tlimit):
			self.myTlimit = args.Tlimit
			self.myTlimit = myTlimit.replace("[","").replace("]","").split(",")
			self.TlimitOpHigh = myTlimit[1].strip().upper()
			self.TlimitOpLow = myTlimit[0].strip().upper()
			allOps = ["LE","LT","GE","GT"]
			if(not(self.TlimitOpLow in allOps) or not(self.TlimitOpHigh in allOps)):
				die("ERROR: on -Tlimit operators must be one of the followings: "+(", ".join(allOps)))
			print "Reading option -Tlimit (Low="+self.TlimitOpLow+", High="+self.TlimitOpHigh+")"
		#determine cooling types
		if(args.cooling):
			myCools = args.cooling.split(",")
			allCools = ["ATOMIC","H2","HD","Z","DH","DUST","H2GP98","COMPTON","CIE",
					"CI","CII","SiI","SiII","OI","OII","FeI","FeII","CONT","CHEM"]
			for coo in myCools:
				if(not(coo in allCools)):
					die("ERROR: Cooling \""+coo+"\" is unknown!\nAvailable coolings are: "+(", ".join(allCools)))

			if("ATOMIC" in myCools): self.useCoolingAtomic = True
			if("H2" in myCools): self.useCoolingH2 = True
			if("H2GP98" in myCools): self.useCoolingH2GP98 = True
			if("HD" in myCools): self.useCoolingHD = True
			if("DH" in myCools): self.useCoolingdH = True
			if("DUST" in myCools): self.useCoolingDust = True
			if("COMPTON" in myCools): self.useCoolingCompton = True
			if("CHEM" in myCools): self.useCoolingChem = True
			if("CIE" in myCools): self.useCoolingCIE = True
			if("CONT" in myCools): self.useCoolingCont = True
			if("Z" in myCools): 
				self.useCoolingZ = self.useCoolingZC = self.useCoolingZCp = self.useCoolingZSi = True
				self.useCoolingZSip = self.useCoolingZO = self.useCoolingZOp = self.useCoolingZFe = True
				self.useCoolingZFep = True
			if("CI" in myCools): self.useCoolingZ = self.useCoolingZC = True
			if("CII" in myCools): self.useCoolingZ = self.useCoolingZCp = True
			if("SiI" in myCools): self.useCoolingZ = self.useCoolingZSi = True
			if("SiII" in myCools): self.useCoolingZ = self.useCoolingZSip = True
			if("OI" in myCools): self.useCoolingZ = self.useCoolingZO = True
			if("OII" in myCools): self.useCoolingZ = self.useCoolingZOp = True
			if("FeI" in myCools): self.useCoolingZ = self.useCoolingZFe = True
			if("FeII" in myCools): self.useCoolingZ = self.useCoolingZFep = True
				
			self.use_cooling = True
			self.hasDust = False
			for aa in argv:
				if("dust=" in aa): self.hasDust = True
			if(self.useCoolingDust and not(self.hasDust)):
				die("ERROR: to include dust cooling you need dust (use -dust=[see help]).")
			if(("CHEM" in myCools) and ("ATOMIC" in myCools)):
				die("ERROR: CHEM and ATOMIC cooling are mutually exclusive!")

			self.use_thermo = True

			print "Reading option -cooling ("+(",".join(myCools))+")"
		
		#determine heating types
		if(args.heating):
			myHeat = args.heating.upper().split(",")
			allHeats = ["COMPRESS","PHOTO","CHEM","DH"]
			for hea in myHeat:
				if(not(hea in allHeats)):
					die("ERROR: Heating \""+hea+"\" is unknown!\nAvailable heatings are: "+(", ".join(allHeats)))

			if("COMPRESS" in myHeat): self.useHeatingCompress = True
			if("PHOTO" in myHeat): self.useHeatingPhoto = True
			if("CHEM" in myHeat): self.useHeatingChem = True
			if("DH" in myHeat): self.useHeatingdH = True

			self.use_thermo = True
			if(not(self.usePhIoniz) and self.useHeatingPhoto):
				die("ERROR: if you use photoheating you have to include potoionization via -usePhIoniz")

			print "Reading option -heating ("+(",".join(myHeat))+")"
	

		#force rwork size
		#if(args.testFile):
		#	myrwork = (arg.strip().replace("-testFile=",""))
		#	self.is_test = True
		#	self.test_name = myrwork
		#	print "Reading option -testFile (file="+str(myrwork)+")"
		#	break

		#force rwork size
		if(args.forceRWORK):
			myrwork = args.forceRWORK
			self.force_rwork = True
			print "Reading option -forceRWORK (RWORK="+str(myrwork)+")"

		#use custom function for coefficient instead of coe_tab()
		if(args.useCustomCoe):
			self.customCoeFunction = args.useCustomCoe
			self.useCustomCoe = True
			print "Reading option -useCustomCoe (Expression="+str(customCoeFunction)+")"

		#use function to append after each ODE
		if(args.useODEConstant):
			self.ODEConstant = args.useODEConstant
			self.useODEConstant = True
			print "Reading option -useODEConstant (Constant="+str(ODEConstant)+")"

		#dust
		hasDustOptions = False
		if(args.dustOptions): hasDustOptions = True
		if(args.dust):
			dustopt = args.dust
			adust = dustopt.split(",")
			self.useDust = True
			if(len(adust)<2): die("ERROR: you must specify dust size and type(s), e.g. -dust=20,C,Si")
			if(self.use_implicit_RHS): die("ERROR: you cannot use dust AND implicit RHS: remove -iRHS option")
			self.dustArraySize = int(adust[0])
			self.dustTypes = adust[1:]
			self.dustTypesSize = len(self.dustTypes)
			print "Reading option -dust (size="+str(self.dustArraySize)+", type(s)="+(",".join(self.dustTypes))+")"
			if(not(hasDustOptions)):
				print "ERROR: -dust flag needs to define -dustOptions=[see help])"
				sys.exit()
		#dust options
		if(args.dustOptions):
			if(not(self.useDust)): die("ERROR: you need -dust=[see help] to activate dust options!")
			dustopt = args.dustOptions
			dustOptions = dustopt.split(",")
			if("GROWTH" in dustOptions): self.useDustGrowth = True
			if("SPUTTER" in dustOptions): self.useDustSputter = True
			if("H2" in dustOptions): self.useDustH2 = True
			if("T" in dustOptions): self.useDustT = True
			print "Reading option -dustOptions (options="+(",".join(dustOptions))+")"

		#project name folder
		if(args.project):
			self.projectName = projectName = args.project
			print "Reading option -project (name="+str(projectName)+")"
			self.buildFolder = "build_"+projectName+"/"
			fout = open(projectName+".kpj","w")
			fout.write((" ".join(argv)))
			fout.close()

		#typeGamma
		if(args.gamma):
			typeGamma = args.gamma
			self.typeGamma = typeGamma.replace("\"","")
			print "Reading option -gamma (gamma="+str(self.typeGamma)+")"

		#ATOL
		if(args.ATOL):
			self.ATOL = args.ATOL
			print "Reading option -atol (atol="+str(self.ATOL)+")"

		#RTOL
		if(args.RTOL):
			self.RTOL = args.RTOL
			print "Reading option -rtol (rtol="+str(self.RTOL)+")"

		#maxord
		if(args.maxord):
			self.maxord = min(max(1,int(args.maxord)),5)
			print "Reading option -maxord (maxord="+str(self.maxord)+")"

		#custom ATOLs
		if(args.customATOL):
			fname = args.customATOL
			print "Reading option -customATOL (file="+fname+")"
			if(not(file_exists(fname))):
				print "ERROR: custom ATOL file \""+fname+"\" does not exist!"
				sys.exit()
			fh = open(fname,"rb")
			for row in fh:
				srow = row.strip()
				if(len(srow)==0): continue
				if(srow[0]=="#"): continue
				arow = [x for x in srow.split(" ") if x.strip()!=""]
				if(len(arow)<2):
					print "ERROR: wrong format in custom ATOL file!"
					print srow
					sys.exit()
				print "ATOL: "+arow[0]+" "+arow[1]
				self.atols.append([arow[0],arow[1]])
			fh.close()



		#custom RTOLs
		if(args.customRTOL):
			fname = args.customRTOL
			print "Reading option -customRTOL (file="+fname+")"
			if(not(file_exists(fname))):
				print "ERROR: custom RTOL file \""+fname+"\" does not exist!"
				sys.exit()
			fh = open(fname,"rb")
			for row in fh:
				srow = row.strip()
				if(len(srow)==0): continue
				if(srow[0]=="#"): continue
				arow = [x for x in srow.split(" ") if x.strip()!=""]
				if(len(arow)<2):
					print "ERROR: wrong format in custom RTOL file!"
					print srow
					sys.exit()
				print "RTOL: "+arow[0]+" "+arow[1]
				self.rtols.append([arow[0],arow[1]])
			fh.close()

		#custom ODEs
		if(args.customODE):
			fname = args.customODE
			print "Reading option -customODE (file="+fname+")"
			if(not(file_exists(fname))):
				print "ERROR: custom ODE file \""+fname+"\" does not exist!"
				sys.exit()
			fh = open(fname,"rb")
			ivarcoe = 0
			for row in fh:
				srow = row.strip()
				if(len(srow)==0): continue
				if(srow[0]=="#"): continue
				#search for variables
				if("@var:" in srow):
					arow = srow.replace("@var:","").split("=")
					if(len(arow)!=2):
						print "ERROR: variable line must be @var:variable=F90_expression"
						print "found: "+srow
						sys.exit()
					print "var: "+arow[0]
					self.coevarsODE[arow[0]] = [ivarcoe,arow[1]]
					ivarcoe += 1 #count variables to sort
					continue #SKIP: a variable line is not a reaction line		
				#search for ODE		
				arow = [x.strip() for x in srow.split("=")]
				if(len(arow)!=2):
					print "ERROR: wrong format in custom ODE file!"
					print srow
					sys.exit()
				print "ODE: "+arow[0]+" "+arow[1]
				self.customODEs.append([arow[0],arow[1]])
			fh.close()


	####################################################
	#load thermochemistry data form chemkin-formatted file
	def load_thermochemistry(self):
		nskip = 99999 #skip comments
		thermo = dict() #prepare dictionary
		fth = open("data/thermo30.dat") #open thermochemistry file
		#loop on file
		for row in fth:
			srow = row.strip()
			if(len(srow)==0): continue #skip empty lines
			arow = srow.split()
			if(arow[0]=="END"): break #break on END
			if(arow[0]=="THERMO"): #start to read data
				nskip = 1 #skip line after thermo
				continue
			#skip comments
			if(nskip>0 or srow[0]=="!"):
				nskip -= 1 #reduce line to be skipped by one
				continue
			#if not number is a species
			if(not(is_number(row[:15]))):
				spec = arow[0].strip().upper() #read species name
				Tmin = arow[len(arow)-4]
				Tmed = arow[len(arow)-2]
				Tmax = arow[len(arow)-3]
				mypoly = [Tmin, Tmed, Tmax] #first data are temperature limits
			else:
				coe = [row[i*15:(i+1)*15] for i in range(5)] #read 5 coefficients
				irow = int(row[5*15:].strip()) #read line number
				mypoly += coe #add coefficients to the coefficients list
				#last line for the given species
				if(irow==4):
					#convert coefficients to floating and skip empty values
					coef = [float(x) for x in mypoly[:17] if x.strip()!=""]
					#check the number of coefficients (3temp+14poly)
					if(len(coef)!=17):
						print "ERROR: wrong format for NASA polynomials!"
						print spec	
						print srow
						sys.exit()
					thermo[spec] = coef #append coefficients to the dictionary
		fth.close()
		self.thermodata = thermo
		print "Thermochemistry data loaded!"



	################################################
	def prepare_massdict(self):
		self.use_RHS_variable = False

		self.separator = "," #separator character

		me =  9.10938188e-28 #electron mass (g)
		mp = 1.67262158e-24 #proton mass (g)
		mn = 1.6725e-24 #neutron mass (g)
		menp = me + mp + mn
		#mass dictionary
		mass_dic = {'H':me+mp,
			'D':menp,
			'He':2.*(menp),
			'Li':3.*(me+mp)+4.*mn,
			'Be':4.*(me+mp)+5.*mn,
			'C':6.*(menp),
			'N':7.*(menp),
			'O':8.*(menp),
			'F':9.*(menp)+mn,
			'Mg':12.*(menp),
			'Na':(me+mp)*11+mn*12,
			'Si':14.*(menp),
			'P':15.*(menp)+mn,
			'S':(menp)*16,
			'Cl':(menp)*17+mn,
			'Fe':(me+mp)*26+mn*29,
			'GRAIN0': 100*6*(menp),
			'GRAIN-': 100*6*(menp)+me,
			'GRAIN+': 100*6*(menp)-me,
			'PAH': 30*6*(menp),
			'PAH-': 30*6*(menp)+me,
			'PAH+': 30*6*(menp)-me,
			'O(1D)':8.*(menp),
			'O(3P)':8.*(menp),
			'CR':0.,
			'M':0.,
			'g':0.,
			'E':me,
			'-':me,
			'+':-me}

		#fake species FK1,FK2,...
		for i in range(10):
			mass_dic['FK'+str(i)] = 1.

		#exited levels of some molecules
		for i in range(9):
			mass_dic['CH2_'+str(i+1)] = 6.*(menp) + 2.*(me+mp)
			mass_dic['SO2_'+str(i+1)] = 16.*(menp) + 2.*8.*(menp)

		#build isotopes (including some non-esistent) as [n]A
		#with -usePlainIsotopes build as nA
		atoms_iso = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Fe","Co","Ni"]
		atoms_p = [i+1 for i in range(20)] + [26,27,28]
		if(len(atoms_iso)!=len(atoms_p)): 
			die("ERROR: in building isotopes the length of the atoms array and the number of protons array mismatch!")
		for aiso in atoms_iso:
			protons = atoms_p[atoms_iso.index(aiso)] #get proton numbers
			for i in range(protons,80):
				iso_name = "["+str(i)+aiso+"]"
				if(self.usePlainIsotopes): iso_name = str(i)+aiso
				mass_dic[iso_name] = protons*(me+mp) + ((i-protons)*mn)

		#sort dictionaries
		self.mass_dic = dict([[k.upper(),v] for (k,v) in mass_dic.iteritems()])
		self.atoms = sorted(mass_dic, key = lambda x: len(x),reverse=True)

	#################################################
	def read_file(self):
		skipDup = self.skipDup
		filename = self.filename
		atoms = self.atoms
		mass_dic = self.mass_dic
		print "Reading from file \""+filename+"\"..."
		spec_names = [] #string
		idx_list = [] #store reaction index in case of -useFileIdx
		pseudo_hash_list = []
		rcount = 0 #count reactions
		found_one = False #flag to control if at least one reaction has been found
		unmatch_idx = False #controls if the reaction index in the file match the sequential index
		#default position in input array line
		iidx = 0 #position of the index
		skipped_dupl = 0 #number of duplicated reaction skipped
		ireact = range(1,4) #positions of the reactants
		iprod = range(4,8) #position of the products
		iTmin = 8 #position of tmin
		iTmax = 9 #position of tmax
		irate = 10 #position of the rate in F90 style
		ivarcoe = 0 #number variable to order dictionary (dictionaries are not ordered by definition!)
		hasFormat = False
		format_items = 4+len(ireact)+len(iprod)
		if(skipDup): fdup = open("duplicates.log","w")
		idxFound = tminFound = tmaxFound = rateFound = True
		specs = []
		reacts = []
		fh = open(filename,"rb") #OPEN FILE
		isComment = False #flag for comment block
		for row in fh:
			srow = row.strip() #stripped row
			if(srow.strip()==""): continue #looks for blank line
			if(srow[0]=="#"): continue #looks for comment line
			if(srow[0:1]=="//"): continue #looks for comment line
			if(srow[0:1]=="/*"): isComment = True #start multiline comment
			#end multiline comment
			if("*/" in srow):
				isComment = False
				continue
			if(isComment): continue #skip if in comment block

			#search for variables
			if("@var:" in srow):
				arow = srow.replace("@var:","").split("=")
				if(len(arow)!=2):
					print "ERROR: variable line must be @var:variable=F90_expression"
					print "found: "+srow
					sys.exit()
				self.coevars[arow[0]] = [ivarcoe,arow[1]]
				ivarcoe += 1 #count variables to sort
				continue #SKIP: a variable line is not a reaction line

			#search for format string
			if("@format:" in srow):
				idxFound = tminFound = tmaxFound = rateFound = False
				hasFormat = True #format flag
				srow = srow.replace("@format:","") #remove 
				print "Found custom format: "+srow
				arow = srow.split(",") #split format line
				#check format (at least 6 elements)
				if(len(arow)<6):
					print "ERROR: format line must contains at least 6 elements"
					print " idx,R,R,P,P,rate"
					print " You provided "+str(len(arow))+" elements:"
					print " "+srow
					sys.exit()
				ipos = 0 #count items in the splitted line
				ireact = [] #reactants index array
				iprod = [] #products index array
				format_items = len(arow)
				#read format elements
				for x in arow:
					x = x.lower().strip() #lower trimmed item
					if(x=="idx"): 
						iidx = ipos #index position
						idxFound = True
					if(x=="r"): ireact.append(ipos) #reactants positions
					if(x=="p"): iprod.append(ipos) #products positions
					if(x=="tmin"):
						iTmin = ipos #min temperature position
						tminFound = True
					if(x=="tmax"): 
						iTmax = ipos #max temperature position
						tmaxFound = True
					if(x=="rate" or x=="k"): 
						irate = ipos #rate in F90 style position
						rateFound = True
					ipos += 1 #increase position
				#check for rate
				if(not(rateFound)):
					print "ERROR: format must contain rate token"
					sys.exit()
			
				continue #SKIP format line (it is not a reaction line)
			arow = srow.split(self.separator,format_items-1) #split only N+1 elements with N seprations
			if(len(arow)!=format_items): 
				print "WARNING: wrong format for reaction "+str(rcount+1)
				print srow
				a = raw_input("Any key to continue q to quit... ")
				if(a=="q"): print sys.exit()
				continue #check line format (N elements, 4=idx+Tmin+Tmax+rate)
			found_one = True
			rcount += 1
			myrea = reaction() #create object reaction

			#use reaction index found into the file
			if(self.useFileIdx):
				if(not(hasFormat) or (hasFormat and idxFound)):
					reaction_idx = int(arow[0]) #index of the reaction (from file)
					if(reaction_idx<1): die("ERROR: reaction index must be > 0!\n Check your reaction file!")
					if(reaction_idx in idx_list):
						die("ERROR: reaction index "+str(reaction_idx)+" already present!\n Check your reaction file!")
					myrea.idx = reaction_idx
					rcount = reaction_idx
			else:
				myrea.idx = rcount
				if(not(hasFormat) or (hasFormat and idxFound)):
					if(rcount!=int(arow[0])): unmatch_idx = True

			reactants = [arow[x] for x in ireact]
			products = [arow[x] for x in iprod]

			opTlist = ["<",">",".LE.",".GE.",".LT.",".GT."]
			myrea.TminOp = self.TlimitOpLow
			myrea.TmaxOp = self.TlimitOpHigh
			for op in opTlist:
				if(not(tminFound)): break
				if(op in arow[iTmin]):
					arow[iTmin] = arow[iTmin].replace(op,"")
					myrea.TminOp = op.replace(">","GT").replace("<","LT").replace(".","")
					break
			for op in opTlist:
				if(not(tmaxFound)): break
				if(op in arow[iTmax]):
					arow[iTmax] = arow[iTmax].replace(op,"")		
					myrea.TmaxOp = op.replace(">","GT").replace("<","LT").replace(".","")
					break

			myrea.Tmin = TminAuto = "2.73d0" #default min temperature
			myrea.Tmax = TmaxAuto = "1.d8" #default max temperature
			if(tminFound): myrea.Tmin = format_double(arow[iTmin]) #get Tmin
			if(tmaxFound): myrea.Tmax = format_double(arow[iTmax]) #get Tmax
			if(tminFound): TminAuto = min(float(arow[iTmin].lower().replace("d","e")), TminAuto)
			if(tmaxFound): TmaxAuto = max(float(arow[iTmax].lower().replace("d","e")), TmaxAuto)

			myrea.krate = arow[irate] #get reaction rate written in F90 style
			if(self.useCustomCoe): myrea.krate = "0.d0" #when custom function is used standard coefficient are set to zero
			#loop over reactants to grep molecules
			for r in reactants:
				if(r.strip()=="g" and not(self.use_photons)): continue
				if(r.strip()!=""):
					mol = parser(r,mass_dic,atoms,self.thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
					myrea.reactants.append(mol) #add molecule object to reactants
			#loop over prodcuts to grep molecules
			for p in products:
				if(p.strip()=="g" and not(self.use_photons)): continue
				if(p.strip()!=""):
					mol = parser(p,mass_dic,atoms,self.thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
					myrea.products.append(mol) #add molecule object to products

			myrea.build_verbatim() #build reaction as string (e.g. A+B->C)
			myrea.reactants = sorted(myrea.reactants, key=lambda r:r.idx) #sort reactants
			myrea.products = sorted(myrea.products, key=lambda p:p.idx) #sort products
			myrea.build_RHS() #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
			myrea.build_phrate() #build photoionization rate
			myrea.check() #check mass and charge conservation
			if(myrea.krate.count("(")!=myrea.krate.count(")")):
				print "ERROR: unbalanced brakets in reaction "+str(myrea.idx)
				print " "+myrea.verbatim
				print " rate = "+myrea.krate
				sys.exit()

			skip_append = False
			if(skipDup):
				myrea.build_pseudo_hash() #build pseudo_hash
				if(myrea.pseudo_hash in pseudo_hash_list): 
					skip_append = True
					skipped_dupl += 1
					fdup.write(str(myrea.idx)+" "+myrea.verbatim+"\n")
					rcount -= 1
				else:
					pseudo_hash_list.append(myrea.pseudo_hash)
			if(not(skip_append)): reacts.append(myrea)

		if(skipDup): 
			fdup.close()
			print "Skipped duplicated reactions:",skipped_dupl

		#check file format
		if(not(found_one)):
			die("ERROR: no valid reactions found in file \""+filename+"\"")
		if(unmatch_idx):
			print "WARNING: index in \""+filename+"\" are not sequential!"

		self.specs = specs
		self.reacts = reacts
		self.TminAuto = TminAuto
		self.TmaxAuto = TmaxAuto

		print "done!"

	###############################################
	def do_reverse(self):
		#do reverse reaction if needed
		reacts = self.reacts
		if(self.useReverse):
			print
			print "reversing reactions..."
			count_reverse = len(reacts)
			reacts_copy = [x for x in reacts] #create a copy of reacts to loop on
			#loop over found reactants
			for myrea in reacts_copy:
				#skip reactions with more than four products
				if(len(myrea.products)>3):
					print "WARNING: in reversing reaction "+myrea.verbatim+" more than 3 products found! Skipped."
					continue
				else:
					count_reverse += 1
					myrev = reaction() #create object reaction for reverse
					myrev.idx = count_reverse
					myrev.Tmin = myrea.Tmin #get Tmin
					myrev.Tmax = myrea.Tmax #get Tmax
					myrev.reactants = myrea.products
					myrev.products = myrea.reactants
					myrev.TminOp = myrea.TminOp #get Tmin operator
					myrev.TmaxOp = myrea.TmaxOp #get Tmax operator
					myrev.build_verbatim() #build reaction as string (e.g. A+B->C)
					myrev.verbatim += " (REV)" #append REV label to reaction string
					myreactants = myrev.reactants #get list of reactants
					myproducts = myrev.products #get list of products
					myrev.krate = myrea.krate #set the forward reaction
					myrev.krate = myrev.doReverse() #compute reverse reaction using Simoncini2013PRL
					myrev.build_RHS() #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
					myrev.check() #check mass and charge conservation
					reacts.append(myrev)
			print "Inverse reaction added: "+str(count_reverse)
			
		self.reacts = reacts
	
	###################################################
	def verifyThermochem(self):
		for x in self.specs:
			if(not(x.name in self.thermodata)):
				print "WARNING: no thermochemical data for "+x.name+"!"



	###########################################add all the metals to cooling
	def addMetals(self):
		Zcools = []
		specs = self.specs
		if(self.useCoolingZC or self.useCoolingZCp): 
			Zcools.append("C") 
			Zcools.append("C+") 
		if(self.useCoolingZO or self.useCoolingZOp): 
			Zcools.append("O") 
			Zcools.append("O+") 
		if(self.useCoolingZSi or self.useCoolingZSip): 
			Zcools.append("Si") 
			Zcools.append("Si+") 
		if(self.useCoolingZFe or self.useCoolingZFep): 
			Zcools.append("Fe")
			Zcools.append("Fe+")
		for zcool in Zcools:
			zFound = False #flag metal found
			#loop on specs to found metals
			for mol in specs:
				#when found break loop
				if(mol.name.lower()==zcool.lower()):
					zFound = True
					break
			#if not found parse and add
			if(not(zFound)):
				print "Adding specie \""+zcool+"\" (request by metal cooling)"
				mymol = parser(zcool,self.mass_dic,self.atoms,self.thermodata)
				mymol.idx = len(specs)+1
				self.specs.append(mymol)
		self.totMetals = "tot_metals = " + (" + ".join(["n(idx_"+x.replace("+","j")+")" for x in Zcools]))
	
	#######################################
	def addReaMin(self):
		for rea in self.reacts:
			if(not("krome_" in rea.krate)): rea.krate = "1d-40 + ("+rea.krate+")"


	##################################
	def computeEnthalpy(self):
		reacts = self.reacts
		for i in range(len(reacts)):
			reacts[i].enthalpy()
			vrb = reacts[i].verbatim.strip()
			avrb = vrb.split("->")
			#print avrb[0]+" "*(15-len(avrb[0])),"->",avrb[1]+" "*(15-len(avrb[1])),reacts[i].dH,"erg"
		self.reacts = reacts
		print "Enthalpy OK!"

	##################################
	def addDust(self):
		dustTypes = self.dustTypes
		specs = self.specs
		dustArraySize = self.dustArraySize
		#add dust to problem
		if(self.useDust):
			print
			#look for dust atoms in species found
			for dType in dustTypes:
				dTypeFound = False #flag atom found
				#loop on specs to found dust atom
				for mol in specs:
					#when found break loop
					if(mol.name.lower()==dType.lower()):
						dTypeFound = True
						break
				#if not found add to specs parsing the name (e.g. Si)
				if(not(dTypeFound)):
						print "Add species \""+dType+"\" (request by dust type)"
						mymol = parser(dType,self.mass_dic,self.atoms,self.thermodata)
						mymol.idx = len(specs)+1
						specs.append(mymol)
			#add dust bins as species
			for dType in dustTypes:
				for i in range(dustArraySize):
						#create the object named dust_type_index
						mymol = molec()
						mymol.name = "dust_"+dType+"_"+str(i)
						mymol.charge = 0
						mymol.mass = 0.e0
						mymol.ename = "dust_"+dType+"_"+str(i)
						mymol.fidx = "idx_dust_"+dType+"_"+str(i)
						mymol.idx = len(specs)+1
						specs.append(mymol)
			print "Dust added:",self.dustArraySize*self.dustTypesSize
		self.specs = specs

	###################################
	def addSpecial(self):
		specs = self.specs
		customODEs = self.customODEs
		#look for photons and CR to add to the species list
		has_g = has_CR = False
		for mol in specs:
			if(mol.name=="g"):
				has_g = True
			if(mol.name=="CR"):
				has_CR = True

		#append custom ODEs as species
		if(len(customODEs)>0):
			for x in customODEs:
				mymol = molec()
				mymol.name = x[0]
				mymol.charge = 0
				mymol.mass = 0.e0
				mymol.ename = mymol.name
				mymol.fidx = "idx_" + mymol.name
				mymol.idx = len(specs)+1
				specs.append(mymol)
				print "Added "+mymol.name+" as species"

		#append CR as species named CR
		if(not(has_CR)):
			mymol = molec()
			mymol.name = "CR"
			mymol.charge = 0
			mymol.mass = 0.e0
			mymol.ename = "CR"
			mymol.fidx = "idx_CR"
			mymol.idx = len(specs)+1
			specs.append(mymol)
			self.photon_species = mymol #store object
			print "Cosmic Rays (CR) added!"

		#append photons as species named g
		if(not(has_g)):
			mymol = molec()
			mymol.name = "g"
			mymol.charge = 0
			mymol.mass = 0.e0
			mymol.ename = "g"
			mymol.fidx = "idx_g"
			mymol.idx = len(specs)+1
			specs.append(mymol)
			self.photon_species = mymol #store object
			print "Photons (g) added!"

		#append temperature as species named Tgas
		mymol = molec()
		mymol.name = "Tgas"
		mymol.charge = 0
		mymol.mass = 0.e0
		mymol.ename = "Tgas"
		mymol.fidx = "idx_Tgas"
		mymol.idx = len(specs)+1
		specs.append(mymol)
		self.Tgas_species = mymol #store object

		#append dummy as species named dummy
		mymol = molec()
		mymol.name = "dummy"
		mymol.charge = 0
		mymol.mass = 0.e0
		mymol.ename = "dummy"
		mymol.fidx = "idx_dummy"
		mymol.idx = len(specs)+1
		specs.append(mymol)
		self.dummy = mymol #store object

		self.specs = specs
	##############################################
	def countSpecies(self):
		#evaluate the number of chemical species (excluding, dust, dummies, Tgas)
		self.nmols = len(self.specs)-4-self.dustArraySize*self.dustTypesSize
		#print found values
		print
		print "ODEs needed:", len(self.specs)
		print "Reactions found:", len(self.reacts)
		print "Species found:", self.nmols

	########################
	def uniq(self,a):
		u = []
		for x in a:
			if(not(x in u)): u.append(x)
		return u

	###############################################
	def dumpNetwork(self):
		#dump species to log file
		fout = open("species.log","w")
		fout.write("#This file contains:\n")
		fout.write("#1. a list of the species used with their indexes\n")
		fout.write("#2. a list of the species to initialize gnuplot\n")
		fout.write("\n")
		idx = 0
		for mol in self.specs:
			idx += 1
			fout.write(str(idx)+"\t"+mol.name+"\t"+mol.fidx+"\n")

		#dump species to gnuplot initialization
		idx = 0
		fout.write("\n")
		fout.write("\n")
		fout.write("#cut-and-paste this into gnuplot to initialize species variables\n")
		fout.write("# nkrome is an optional offset\n")
		fout.write("nkrome = 0\n")
		inits = []
		for mol in self.specs:
			idx += 1
			inits.append("krome_"+mol.fidx+" = "+str(idx)+" + nkrome")
		fout.write(("\n".join(inits))+"\n")
		fout.close()	
	
		#dump reactions to log file
		fout = open("reactions.log","w")
		idx = 0
		for rea in self.reacts:
			idx += 1
			fout.write(str(rea.idx)+"\t"+rea.verbatim+"\n")
		fout.close()	
		
		#dump network to dot file
		fout = open("network.dot","w")
		ntw = dict()
		dot = "digraph{\n"
		for rea in self.reacts:
			react = self.uniq(rea.reactants)		
			prods = self.uniq(rea.products)
			for x in react:
				dot += "\""+x.name+"\" -> k"+str(rea.idx) +"\n"
			for y in prods:
				dot += "k"+str(rea.idx) +" -> \""+y.name+"\";\n"
		dot +="}\n"
		fout.write(dot)
		fout.close
	
	##############################################
	def simpleCHeader(self):
		if(not(self.wrapC)): return
		fout = open("krome.h","w")
		fout.write("#ifndef KROME_H_\n")
		fout.write("#define KROME_H_\n")
		fout.write("\n")
		idx = 0
		for mol in self.specs:
			idx += 1
			fout.write("extern int krome_"+mol.fidx+" = "+str(idx-1)+"; //"+mol.name+"\n")
		fout.write("\n")
		fout.write("extern char* krome_names[] = {\n")
		fout.write(",\n".join(["  \""+x.name+"\"" for x in self.specs])+"\n};")
		fout.write("\n")
		fout.write("extern int krome_nmols = "+str(self.nmols)+"; //number of species\n")
		fout.write("extern int krome_nrea = "+str(len(self.reacts))+"; //number of reactions\n")
		fout.write("\n")
		fout.write("#endif\n")
		fout.close()

	###############################################
	def showODE(self):
		dustArraySize = self.dustArraySize
		specs = self.specs
		nmols = self.nmols
		#include dust in ODE array partition
		sdust = sdustT = ""
		if(self.useDust):
			#include dust by types
			for dType in self.dustTypes:
				sdust += " + ["+str(dustArraySize)+" "+dType+"-dust]"
		sODEpart = "ODE partition: [" + str(nmols) + " atom/mols] "+sdust+sdustT+" + [1 CR] + [1 PHOT] + [1 Tgas] + [1 dummy] = "
		sODEpart += str(len(specs))+" ODEs"
		print sODEpart

		#if species are few print list of species
		if(len(specs)<20): print "ODEs list: "+(", ".join([x.name for x in specs]))
		print "Species list in species.log"
		print

	############################################
	def createODE(self):
		reacts = self.reacts
		specs = self.specs
		dustTypes = self.dustTypes
		dustArraySize = self.dustArraySize
		dustTypesSize = len(dustTypes)
		ndust = dustArraySize*dustTypesSize
		nmols = self.nmols
		dummy = self.dummy
		#create explicit differentials
		dns = ["dn("+str(sp.idx)+") = 0.d0" for sp in specs] #initialize
		for rea in reacts:
			for r in rea.reactants:
				dns[r.idx-1] = dns[r.idx-1].replace(" = 0.d0"," =")
				dns[r.idx-1] += " -"+rea.RHS
			for p in rea.products:
				dns[p.idx-1] = dns[p.idx-1].replace(" = 0.d0"," =")
				dns[p.idx-1] += " +"+rea.RHS

		#add dust to ODEs
		if(self.useDust):
			j = 0
			for dType in dustTypes:
				for i in range(dustArraySize):
					myndust = dustArraySize*dustTypesSize
					j += 1
					#******growth / sputtering******
					if(self.useDustGrowth):
						dns[nmols+j-1] += " + krome_dust_grow(n("+str(nmols+j)+"),n(idx_"+dType+"),Tgas"
						dns[nmols+j-1] += ",krome_dust_T("+str(j)+"),vgas,krome_dust_asize("+str(j)+"))"
					if(self.useDustSputter):
						dns[nmols+j-1] += " &\n- krome_dust_sput(Tgas,krome_dust_asize("
						dns[nmols+j-1] += str(j)+"),ntot,n("+str(nmols+j)+"))"

		#find the maximum number of products and reactants
		maxnprod = maxnreag = 0
		for rea in reacts:
			maxnprod = max(maxnprod, len(rea.products))
			maxnreag = max(maxnreag, len(rea.reactants))

		print "Max number of reactants:", maxnreag
		print "Max number of products:", maxnprod

		#create implicit RHS arrays
		arr_rr = [[] for i in range(maxnreag)]
		arr_pp = [[] for i in range(maxnprod)]
		for rea in reacts:
			for i in range(maxnreag):
				if(i<len(rea.reactants)): arr_rr[i].append(rea.reactants[i].idx)
				else: arr_rr[i].append(dummy.idx)
			for i in range(maxnprod):
				if(i<len(rea.products)): arr_pp[i].append(rea.products[i].idx)
				else: arr_pp[i].append(dummy.idx)

		self.maxnreag = maxnreag
		self.maxnprod = maxnprod
		implicit_arrays = ""
		for i in range(len(arr_rr)):
			implicit_arrays += "arr_r"+str(i+1)+"(:) = (/"+(",".join([str(x) for x in arr_rr[i]]))+"/)\n"
		for i in range(len(arr_pp)):
			implicit_arrays += "arr_p"+str(i+1)+"(:) = (/"+(",".join([str(x) for x in arr_pp[i]]))+"/)\n"
		self.implicit_arrays = implicit_arrays

		#wrap RHS (e.g. knn+knn=2*knn)
		dnw = []
		idn = 0
		ldns = dns
		for dn in ldns:
			if(idn>nmols):
				dnw.append(dn)
				continue
			RHSs = []
			RHSc = []
			dnspl =  dn.strip().split()
			dns = " ".join(dnspl[:2])
			for p in dnspl[2:]:
				if(not(p in RHSs)):
					RHSs.append(p)
					RHSc.append(0)
				RHSc[RHSs.index(p)] += 1
			for i in range(len(RHSs)):
				#if((i+1) % 4 == 0): dns += "&\n" #break long lines
				if("-" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("-"," &\n-"+str(RHSc[i])+".d0*")
				if("+" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("+"," &\n+"+str(RHSc[i])+".d0*")
				if(RHSc[i]==1): dns+= " &\n"+RHSs[i]
			#dns = dns.replace("*"," * ")
			dnw.append(dns)
			idn += 1
		#dnw = [x.replace("+"," &\n+").replace("-"," &\n-") for x in ldns]
		self.dnw = dnw


	###############################################
	#build the jacobian
	def createJAC(self):
		reacts = self.reacts
		specs = self.specs
		Tgas_species = self.Tgas_species
		ndust = self.dustArraySize * self.dustTypesSize
		nmols = self.nmols
		#create explicit JAC for chemical species
		neq = len(specs)
		jac = [["pdj("+str(j+1)+") = 0.d0" for i in range(neq)] for j in range(neq)] #init jacobian
		jsparse = [[0 for i in range(neq)] for j in range(neq)] #sparsity matrix
		for rea in reacts:
			for ri in range(len(rea.reactants)):
				r1 = rea.reactants[ri]
				sjac = "k("+str(rea.idx)+")"
				for rj in range(len(rea.reactants)):
					r2 = rea.reactants[rj]
					if(ri!=rj):
						sjac += "*n("+str(r2.fidx)+")"
				for rr in rea.reactants:
					jac[rr.idx-1][r1.idx-1] = jac[rr.idx-1][r1.idx-1].replace(" = 0.d0"," =")
					jac[rr.idx-1][r1.idx-1] += " &\n-"+sjac
					jsparse[rr.idx-1][r1.idx-1] = 1
				for pp in rea.products:
					jac[pp.idx-1][r1.idx-1] = jac[pp.idx-1][r1.idx-1].replace(" = 0.d0"," =")
					jac[pp.idx-1][r1.idx-1] += " &\n+"+sjac
					jsparse[pp.idx-1][r1.idx-1] = 1

		#this part compress jacobian terms (e.g. k1*H*H+k1*H*H => 2*k1*H*H)
		jacnew = [["" for i in range(neq)] for j in range(neq)] #create empty jacobian neq*neq
		#loop over jacobian
		for i in range(neq):
			for j in range(neq):
				jcount = dict() #dictionary to count items
				jpz = jac[i][j] #jacobian element
				ajpz = [x.strip() for x in jpz.split("&\n")] #split on returns
				#loop over pieces and count (first piece is initialization)
				for pz in ajpz[1:]:
					if(pz in jcount):
						jcount[pz] += 1 #already found, add 1
					else:
						jcount[pz] = 1 #newly found, init to 1
				jacel = ajpz[0] #element start
				#loop on counted pieces
				for (k,v) in jcount.iteritems():
					if(v>1):
						jacel += " "+k[0]+str(v)+".d0*"+k[1:] #add multiplication factor
					else:
						jacel += " "+k #add piece
				#print jacel
				jacel = jacel.replace("+"," &\n +").replace("-"," &\n -") #new line for each term
				jacnew[i][j] = jacel #update element into the jacobian
		jac = jacnew #replace old jacobian

		#create approximated Jacobian term for Tgas
		for i in range(neq):
			if(not(self.use_thermo)): break
			jsparse[i][Tgas_species.idx-1] = jsparse[Tgas_species.idx-1][i] = 1
			s = str(i+1)
			if(specs[i].name!="dummy" and specs[i].name!="CR" and specs[i].name!="g"):
				jac[Tgas_species.idx-1][i] = "dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &\n"
				jac[Tgas_species.idx-1][i] += "* (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))\n"
				jac[Tgas_species.idx-1][i] += "nn(:) = n(:)\n"
				if(self.deltajacMode=="RELATIVE"):
					jac[Tgas_species.idx-1][i] += "dnn = n("+s+")*"+str(self.deltajac)+"\n"
				elif(self.deltajacMode=="ABSOLUTE"):
					jac[Tgas_species.idx-1][i] += "dnn = "+str(self.deltajac)+"\n"
				else:
					die("ERROR: deltajacMode unknonw!"+self.deltajacMode)
				jac[Tgas_species.idx-1][i] += "if(dnn>0.d0) then\n"
				jac[Tgas_species.idx-1][i] += " nn("+s+") = n("+s+") + dnn\n"
				jac[Tgas_species.idx-1][i] += " dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &\n"
				jac[Tgas_species.idx-1][i] += "  * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))\n"
				jac[Tgas_species.idx-1][i] += " pdj(idx_Tgas) = (dn1-dn0)/dnn\n"
				jac[Tgas_species.idx-1][i] += "end if\n"
				#jac[Tgas_species.idx-1][i] = "if(abs(n("+s+") - jac_nold(" + s 
				#jac[Tgas_species.idx-1][i] += "))>1d-10) pdj(idx_Tgas) = (jac_dn(idx_Tgas) - jac_dnold(idx_Tgas)) / (n("
				#jac[Tgas_species.idx-1][i] += s + ") - jac_nold(" + s + "))"

		self.jac = jac
		self.jsparse = jsparse

	#######################################
	def IACJAC(self):
		neq = len(self.specs)
		jsparse = self.jsparse
		#create IAC/JAC (see DLSODES manual)
		isum = 1
		ia = [1]
		ja = []
		for i in range(neq):
			isum += sum(jsparse[i])
			ia.append(isum)
			for j in range(neq):
				if(jsparse[i][j]==1): ja.append(j+1)
		nnz = len(ja)
		iaf = "IWORK(31:" + str(30+neq+1) + ") = (/" + (",".join([str(x) for x in ia])) + "/)"
		jaf = "IWORK("+str(31+neq+1)+":"+str(31+neq+nnz)+") = (/"+(",".join([str(x) for x in ja]))+"/)"
		pnnz = round(nnz*100./neq/neq,2)
		print
		print "Jacobian non-zero elements:",nnz,"over",neq*neq
		print "("+str(pnnz)+"% of total elements, sparsity = "+str(100.-pnnz)+"%)"

		#if MF=222 no need for sparsity structure arrays
		if(self.solver_MF == 222):
			iaf = jaf = ""
		self.ja = ja
		self.ia = ia
		self.iaf = iaf
		self.jaf = jaf
	####################################
	def solverParams(self):
		specs = self.specs
		reacts = self.reacts

		solver_MF = self.solver_MF
		solver_moss = int(solver_MF/100)
		solver_meth = int(solver_MF/10) % 10
		solver_miter = solver_MF % 10

		print "solver info:"
		print " MF:",solver_MF
		print " MOSS+METH+MITER:","+".join([str(x) for x in [solver_moss,solver_meth,solver_miter]])


		#estimate size of RWORK array (see DLSODES manual)
		neq = len(specs) #number of equations
		nrea = len(reacts) #number of reactions
		lenrat = 2
		nnz_estimate = pow(neq,2) - nrea #estimate non-zero elements (inaccurate)
		nnz = len(self.ja)
		nnz = max(nnz, nnz_estimate) #estimate may be more realistic than actual nnz calculation!
		if(nnz<0):
			die("ERROR: nnz<0, please check RWORK size estimation!")

		#calculate LWM
		if(solver_miter==0):
			lwm = 0
		elif(solver_miter==1):
			lwm = 2*nnz + 2*neq + (nnz+9*neq)/lenrat
		elif(solver_miter==2):
			lwm = 2*nnz + 2*neq + (nnz+10*neq)/lenrat
		elif(solver_miter==3):
			lwm = neq + 2
		else:
			die("ERROR: solver_miter value "+str(solver_miter)+" unknown!")

		#calculate LRW as in DLSODES documentation
		if(solver_MF==10):
			lrw = 20+16*neq
		elif(solver_MF in [11,111,211,12,112,212]):
			lrw = 20+16*neq+lwm
		elif(solver_MF==13):
			lrw = 20+17*neq
		elif(solver_MF==20):
			lrw = 20+9*neq
		elif(solver_MF in [21,121,221,22,122,222]):
			lrw = 20+9*neq+lwm
		elif(solver_MF==23):
			lrw = 22+10*neq
		else:
			die("ERROR: solver_MF value "+str(solver_MF)+" unknown in LRW calculation!")

		if(self.force_rwork):
			lrw = myrwork
		#lrw = int(20+(2+0.5)*nnz + (11+4.5)*neq) #RWORK size

		print " LWM:",lwm,"LRW:",lrw
		self.lwm = lwm
		self.lrw = lrw

	#######################################
	def prepareBuild(self):
		buildFolder = self.buildFolder
		#create build folder if not exists
		if(not(os.path.exists(buildFolder))):
			os.mkdir(buildFolder)
			print "Created "+buildFolder
		#remove everything
		if(self.cleanBuild):
			clear_dir(buildFolder) #clear the build directory
			print "Deleted all contents in "+buildFolder
		#remove only selected krome files
		else:
			underFiles = ["commons","constants","cooling","dust","ode","reduction","subs","tabs","user",]
			delFiles = [buildFolder+"krome_"+x+".f90" for x in underFiles] + [buildFolder+x for x in ["krome.f90","opkda1.f","opkda2.f"]]
			delFiles += glob.glob(buildFolder+"*~") + glob.glob(buildFolder+"*.mod")
			delFiles += glob.glob(buildFolder+"*_genmod.f90") + glob.glob(buildFolder+"*.i90")
			#for dfile in delFiles:
				#print "deleting "+dfile


		print
		print "Prepearing files in /build..."

		#build all the modules in a single file krome_all.f90
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","w")
			fout.close()

	################################################
	def makeCommons(self):

		reacts = self.reacts
		specs = self.specs
		buildFolder = self.buildFolder
		#*********COMMONS****************
		#write parameters in krome_commons.f90
		print "- writing krome_commons.f90...",
		fh = open("src/krome_commons.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_commons.f90","w")

		#photoionization variables and functions
		phvars = []
		pheatvars = []
		if(self.usePhIoniz):
			for react in reacts:
				if("krome_kph_" in react.krate):
					phvars.append(react.krate)
					pheatvars.append(react.krate.replace("krome_kph_","krome_pheat_"))

		#common dust optical variables
		optVariables = ""
		if(self.useDust):
			for dType in self.dustTypes:
				optVariables += "real*8,allocatable::dust_opt_asize_"+dType+"(:), dust_opt_nu_"
				optVariables += dType+"(:),dust_opt_Qabs_"+dType+"(:,:)\n"
				optVariables += "real*8,allocatable::dust_opt_Em_"+dType+"(:,:),dust_opt_Tbb_"+dType+"(:)\n"

		#common variables
		skip = False
		for row in fh:
			srow = row.strip()
			if(srow == "#KROME_species_index"):
				for x in specs:
					fout.write("\tinteger,parameter::" + x.fidx + "=" + str(x.idx) + "\n")
			elif(srow == "#KROME_parameters"):
					ndust = self.dustArraySize*self.dustTypesSize
					fout.write("\tinteger,parameter::nrea=" + str(len(self.reacts)) + "\n")
					fout.write("\tinteger,parameter::nmols=" + str(self.nmols) + "\n")
					fout.write("\tinteger,parameter::nspec=" + str(len(specs)) + "\n")
					fout.write("\tinteger,parameter::ndust=" + str(ndust) + "\n")
					fout.write("\tinteger,parameter::ndustTypes=" + str(self.dustTypesSize) + "\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			elif(srow == "#KROME_implicit_arr_r"):
				for j in range(self.maxnreag):
					fout.write("integer::arr_r"+str(j+1)+"(nrea)\n")
			elif(srow == "#KROME_implicit_arr_p"):
				for j in range(self.maxnprod):
					fout.write("integer::arr_p"+str(j+1)+"(nrea)\n")
			elif(srow == "#KROME_photo_variables" and self.usePhIoniz and len(phvars)>0):
				fout.write("real*8::"+(",".join(phvars))+"\n")
			elif(srow == "#KROME_photoheating_variables" and self.useHeatingPhoto):
				fout.write("real*8::"+(",".join(pheatvars))+"\n")
			elif(srow == "#KROME_opt_variables"):
				fout.write(optVariables)
			else:
				if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"


	##################################
	def makeConstants(self):
		buildFolder = self.buildFolder
		constants = []
		constants.append(["boltzmann_eV", "8.617332478d-5","eV / K"]) 
		constants.append(["boltzmann_J", "1.380648d-23","J / K"])
		constants.append(["boltzmann_erg", "1.380648d-16","erg / K"]) 
		constants.append(["planck_eV","4.135667516d-15","eV s"]) 
		constants.append(["planck_J","6.62606957d-34","J s"])  
		constants.append(["planck_erg","6.62606957d-27","erg s"]) 
		constants.append(["gravity","6.674d-8","cm3 / g / s2"])      
		constants.append(["e_mass","9.10938188d-28","g"]) 
		constants.append(["p_mass","1.67262158d-24","g"]) 
		constants.append(["clight","2.99792458e10","cm/s"]) 
		constants.append(["pi","3.14159265359d0","#"]) 
		constants.append(["eV_to_erg","1.60217646d-12","eV -> erg"]) 
		constants.append(["seconds_per_year","365d0*24d0*3600d0","yr -> s"]) 
		constants.append(["kvgas_erg","8.d0*boltzmann_erg/pi/p_mass",""]) 
		constants.append(["pre_planck","2.d0*planck_erg/clight**2","erg/cm2*s3"]) 
		constants.append(["exp_planck","planck_erg / boltzmann_erg","s*K"]) 
		constants.append(["stefboltz_erg","5.670373d-5","erg/s/cm2/K4"])
		constants.append(["N_avogadro","6.0221d23","#"]) 
		constants.append(["Rgas_J","8.3144621d0","J/K/mol"]) 
		constants.append(["Rgas_kJ","8.3144621d-3","kJ/K/mol"])



		#********* CONSTANTS ****************
		fh = open("src/krome_constants.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_constants.f90","w")
		
		#prepares list of constants
		const = "!constants\n"
		for x in constants:
			const += "real*8,parameter::" + x[0] + " = " + x[1] + " !" + x[2] + "\n"

		#replace pragmas
		for row in fh:
			if(row.strip()=="#KROME_constant_list"):
				fout.write(const)
			if(row[0]!="#"):
				fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		self.constantList = constants


	########################################
	def makeUserCommons(self):
		buildFolder = self.buildFolder
		#*********USER COMMONS****************
		#write parameters in krome_user_commons.f90
		if(not(file_exists(buildFolder+"krome_user_commons.f90"))):
			print "- writing krome_user_commons.f90...",

			fh = open("src/krome_user_commons.f90")
			fouta = open(self.buildFolder+"krome_user_commons.f90","w")

			for row in fh:
				row = row.replace("#KROME_header", get_licence_header(self.version, self.codename))

				if(row[0]!="#"): fouta.write(row)

			fouta.close()
			print "done!"
		else:
			print "WARNING: krome_user_commons.f90 already found in "+buildFolder+" : not replaced!"

	###################################################
	def makeSubs(self):
		buildFolder = self.buildFolder
		reacts = self.reacts
		specs = self.specs
		thermodata = self.thermodata
		coevars = self.coevars
		#*********SUBS****************
		#write parameters in krome_subs.f90
		print "- writing krome_subs.f90...",
		fh = open("src/krome_subs.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_subs.f90","w")

		#create list of temperature shortcuts
		sclist = []
		for rea in reacts:
			sclist = get_Tshortcut(rea,sclist)

		#loop on src file and replace pragmas
		for row in fh:
			srow = row.strip()
			if(srow == "#KROME_krates"):
				for x in reacts:
					sTlimit = ""
					if(x.kphrate==None and self.useTlimits): 
						sTlimit = "if(Tgas."+x.TminOp+"."+x.Tmin
						sTlimit += " .and. Tgas."+x.TmaxOp+"."+x.Tmax+")"
					kstr = "!" + x.verbatim+"\n"
					kstr += "\t" + sTlimit + " k("+str(x.idx)+") = " + x.krate
					kstr = truncF90(kstr, 60,"*")
					fout.write(truncF90(kstr, 60,"/")+"\n\n")
			elif(srow == "#KROME_implicit_arrays"):
				fout.write(truncF90(self.implicit_arrays,60,","))
			elif(srow == "#KROME_initcoevars"):
				if(len(coevars)==0): continue
				kvars = "real*8::"+(",".join([x for x in coevars.keys()]))
				fout.write(kvars+"\n")
			elif(srow == "#KROME_coevars"):
				if(len(coevars)==0): continue
				klist = [[k+" = "+v[1]+"\n",v[0]] for k,v in coevars.iteritems()] #this mess is to sort dict
				klist = sorted(klist, key=lambda x: x[1])
				fout.write("".join([x[0] for x in klist]))
			elif(srow == "#KROME_masses"):
				for x in specs:
					massrow = "\tget_mass("+str(x.idx)+") = " + str(x.mass).replace("e","d") + "\t!" + x.name + "\n"
					fout.write(massrow.replace("0.0","0.d0"))
			elif(srow == "#KROME_names"):
				for x in specs:
					fout.write("\tget_names("+str(x.idx)+") = \"" + x.name + "\"\n")
			elif(srow == "#KROME_charges"):
				for x in specs:
					fout.write("\tget_charges("+str(x.idx)+") = " + str(x.charge) + ".d0 \t!" + x.name + "\n")
			elif(srow == "#KROME_reaction_names"):
				for x in reacts:
					kstr = "\tget_rnames("+str(x.idx)+") = \"" + x.verbatim +"\""
					fout.write(kstr+"\n")
			elif(srow == "#KROME_Tshortcuts"):
				for shortcut in sclist:
					fout.write(shortcut+"\n")
			elif(srow == "#KROME_rvars"):
				if(self.maxnreag>0):
					fout.write("integer::"+(",".join(["r"+str(j+1) for j in range(self.maxnreag)]))+"\n")
			elif(srow == "#KROME_arrs"):
				if(self.maxnreag>0):
					for j in range(self.maxnreag):
						fout.write("r"+str(j+1)+" = arr_r"+str(j+1)+"(i)\n")
			elif(srow == "#KROME_arr_flux"):
				if(self.maxnreag>0):
					fout.write("arr_flux(i) = k(i)*"+("*".join(["n(r"+str(j+1)+")" for j in range(self.maxnreag)]))+"\n")
			elif(srow == "#KROME_sum_H_nuclei"):
				hsum = []
				for x in specs:
					if(not("H") in x.atomcount): continue
					if(x.atomcount["H"]==0): continue
					hmult = ("*"+format_double(x.atomcount["H"]) if x.atomcount["H"]>1 else "")
					hsum.append("n("+x.fidx+")"+hmult)
				if(len(hsum)==0): hsum.append("0.d0")
				fout.write("nH = "+(" + &\n".join(hsum))+"\n")
			elif(srow == "#KROME_var_reverse"):
				slen = str(len(specs))
				fout.write("real*8::p1("+slen+",7), p2("+slen+",7), Tlim("+slen+",3), p(7)\n")
			elif(srow == "#KROME_kc_reverse"):
				datarev = ""
				sp1 = sp2 = spt = ""
				for x in specs:
					if(min(x.poly1)==0 and max(x.poly1)==0): continue
					sp1 += "p1("+x.fidx+",:)  = (/" + (",&\n".join([format_double(pp) for pp in x.poly1])) + "/)\n"
					sp2 += "p2("+x.fidx+",:)  = (/" + (",&\n".join([format_double(pp) for pp in x.poly2])) + "/)\n"
					spt += "Tlim("+x.fidx+",:)  = (/" + (",&\n".join([format_double(pp) for pp in x.Tpoly])) + "/)\n"
				fout.write(sp1+sp2+spt)
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			else:
				fout.write(row)
		if(not(self.buildCompact)):
			fout.close()
		print "done!"


	##############################################
	def makePhoto(self):

		buildFolder = self.buildFolder
		reacts = self.reacts
		#********* PHOTO ****************
		print "- writing krome_photo.f90...",
		fh = open("src/krome_photo.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_photo.f90","w")

		phvars = []
		pheatvars= []
		ph_func = ph_qromos = ph_heat = ph_heat_qromos = ph_zero = ph_heat_print = ""
		for react in reacts:
			phstuff = get_ph_stuff(react)
			if(phstuff==None): continue
			ph_func += phstuff["ph_func"]+"\n"
			ph_qromos += phstuff["qromos"]+"\n"
			ph_heat += phstuff["ph_heat"]+"\n"
			ph_name = phstuff["reaname"]

			phvars.append("krome_kph_"+ph_name)
			pheatvars.append("krome_pheat_"+ph_name)

			#string to print computed photoheating values (uncomment for debug)
			#ph_heat_print += "print '(a20,E11.3,a1,a6)',\"" + react.verbatim + "\", " + "krome_pheat_"+ ph_name + ",\"\",\"eV/s\"\n"

			#string containing photoheating computations (integrals)
			ph_heat_qromos += "\t" + phstuff["qromos"].replace("sigma_", "heat_").replace("krome_kph_","krome_pheat_") + "\n"

		#inizialization of photoionization and photoheating common variables to zero
		for x in phvars:
			ph_zero += x + " = 0.d0\n"
		for x in pheatvars:
			ph_zero += x + " = 0.d0\n"

		#replace photoionization and photoheating functions
		skip = False
		for row in fh:
			if(row.strip() == "#IFKROME_usePhIoniz" and not(self.usePhIoniz)): skip = True
			if(row.strip() == "#ENDIFKROME"): skip = False

			if(skip): continue
			#replace krome variables
			if(self.usePhIoniz):
				row = row.replace("#KROME_photo_functions", ph_func +"\n")
				row = row.replace("#KROME_photo_qromos", ph_qromos +"\n")
				row = row.replace("#KROME_photo_init_zero", ph_zero +"\n")
			if(self.useHeatingPhoto):
				row = row.replace("#KROME_photo_heating_qromos", ph_heat_qromos +"\n")
				row = row.replace("#KROME_photo_heating_functions", ph_heat +"\n")
				row = row.replace("#KROME_photo_heating_print", ph_heat_print +"\n")


			if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"

	########################################################
	def makeTabs(self):
		buildFolder = self.buildFolder
		#********* TABS ****************
		print "- writing krome_tabs.f90...",
		fh = open("src/krome_tabs.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_tabs.f90","w")


		skip = False
		for row in fh:
			if(row.strip() == "#IFKROME_useCustomCoe" and not(self.useCustomCoe)): skip = True
			if(row.strip() == "#ENDIFKROME"): skip = False

			if(row.strip() == "#IFKROME_useTabs" and not(self.useTabs)): skip = True
			if(row.strip() == "#ENDIFKROME"): skip = False

			if(row.strip() == "#IFKROME_useStandardCoe" and ((self.useCustomCoe) or (self.useTabs))): skip = True
			if(row.strip() == "#ENDIFKROME"): skip = False

			if(skip): continue
			row = row.replace("#KROME_logTlow", "ktab_logTlow = log10(max("+str(self.TminAuto)+",2.73d0))")
			row = row.replace("#KROME_logTup", "ktab_logTup = log10("+str(self.TmaxAuto)+")")
			if(self.useCustomCoe): row = row.replace("#KROMEREPLACE_customCoeFunction", self.customCoeFunction)

			if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"




	###############################
	def makeDust(self):
		buildFolder = self.buildFolder
		useDustT = self.useDustT
		#********* DUST ****************
		print "- writing krome_dust.f90...",
		fh = open("src/krome_dust.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_dust.f90","w")


		dustPartnerIdx = dustQabs = dustOptInt = "" 
		if(self.useDust):
			itype = 0 #dust type index
			#loop in dust types
			for dType in self.dustTypes:
				itype += 1 #increase index
				dustPartnerIdx += "krome_dust_partner_idx("+str(itype)+") = idx_"+dType+"\n"
				if(useDustT):
					dustQabs += "call init_Qabs(\"opt"+dType+".dat\",dust_opt_Qabs_"+dType
					dustQabs += ",dust_opt_asize_"+dType+", &\n dust_opt_nu_"+dType+")\n"
					dustOptInt += "call dustOptIntegral(dust_opt_Em_"+dType+",dust_opt_Tbb_"+dType+","
					dustOptInt += "dust_opt_asize_"+dType+",&\n dust_opt_nu_"+dType+", dust_opt_Qabs_"+dType+")\n"


		skip = False
		for row in fh:
			srow = row.strip()
			if(srow == "#IFKROME_useDust" and not(self.useDust)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			row = row.replace("#KROME_dustPartnerIndex", dustPartnerIdx)
			row = row.replace("#KROME_init_Qabs", dustQabs)
			row = row.replace("#KROME_opt_integral", dustOptInt)
		 
			if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"


	################################
	def makeCooling(self):
		buildFolder = self.buildFolder
		reacts = self.reacts
		specs = self.specs
		#*********COOLING****************
		#write header in krome_cooling.f90
		print "- writing krome_cooling.f90...",
		fh = open("src/krome_cooling.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_cooling.f90","w")


		#create coefficients, flux, and cooling functions for dH cooling
		i = 0
		dH_varsa = []
		dH_coe = dH_cool = ""
		if(self.useCoolingdH):
			for rea in reacts:
				if(rea.dH!=None and rea.dH<0.e0):
					i += 1 #count cooling reactions
					kvar = "k"+str(i) #local variable for coefficient
					dH_varsa.append(kvar) #variables array
					dH_coe += "!"+rea.verbatim+"\n" #reaction comment
					dH_coe += kvar+" = 0.d0\n"
					dH_coe += "if(Tgas."+self.TlimitOpLow+"."+rea.Tmin+" .and. Tgas."+self.TlimitOpHigh+"."+rea.Tmax+") then\n"
					dH_coe += kvar+" = "+rea.krate+"\n" #evaluate reation coefficient
					dH_coe += "end if\n\n" #evaluate reation coefficient
					dH_cool += "cool = cool + "+kvar+"*"+("*".join(["n("+x.fidx+")" for x in rea.reactants])) #evautate cooling
					dH_cool += "*"+(str(abs(rea.dH))).replace("e","d")+"\n"


		#bremsstrahlung for all the ions as charge**2*n_ion
		skip = False
		bms_ions = "bms_ions ="
		for x in specs:
			charge = x.charge
			if(charge>0):
				mult = ""
				if(charge>1): mult = str(charge*charge)+".d0*"
				bms_ions += " +"+mult+"n("+x.fidx+")"

		useCoolingZ = self.useCoolingZ
		#loop on source to replace pragmas
		for row in fh:
			if(row.strip() == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			else:
				srow = row.strip()
				#enthalpic
				if(srow == "#IFKROME_useCoolingdH" and (not(self.useCoolingdH) or len(dH_varsa)==0)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#dust
				if(srow == "#IFKROME_useCoolingDust" and not(self.useCoolingDust)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#metals
				if(srow == "#IFKROME_useCoolingZ" and not(useCoolingZ)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZ"): skip = False

				#individual metals
				if(srow == "#IFKROME_useCoolingZC" and not(self.useCoolingZC)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZC" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZCp" and not(self.useCoolingZCp)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZCp" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZSi" and not(self.useCoolingZSi)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZSi" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZSip" and not(self.useCoolingZSip)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZSip" and useCoolingZ): skip = False

				if(srow == "#IFKROME_useCoolingZO" and not(self.useCoolingZO)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZO" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZOp" and not(self.useCoolingZOp)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZOp" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZFe" and not(self.useCoolingZFe)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZFe" and useCoolingZ): skip = False
		
				if(srow == "#IFKROME_useCoolingZFep" and not(self.useCoolingZFep)): skip = True
				if(srow == "#ENDIFKROME_useCoolingZFep" and useCoolingZ): skip = False


				#CEN
				if(srow == "#IFKROME_useCoolingAtomic" and not(self.useCoolingAtomic)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#H2
				if(srow == "#IFKROME_useCoolingH2" and not(self.useCoolingH2)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#H2GP92
				if(srow == "#IFKROME_useCoolingH2GP" and not(self.useCoolingH2GP98)): skip = True
				if(srow == "#ENDIFKROME"): skip = False
		
				#HD
				if(srow == "#IFKROME_useCoolingHD" and not(self.useCoolingHD)): skip = True
				if(srow == "#ENDIFKROME"): skip = False
		
				#COMPTON
				if(srow == "#IFKROME_useCoolingCompton" and not(self.useCoolingCompton)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#CIE
				if(srow == "#IFKROME_useCoolingCIE" and not(self.useCoolingCIE)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				#Continuum
				if(srow == "#IFKROME_useCoolingContinuum" and not(self.useCoolingCont)): skip = True
				if(srow == "#ENDIFKROME"): skip = False

				if(skip): continue
		
				#replace pragma for total metals
				row = row.replace("#KROME_tot_metals", self.totMetals)
				
				if(self.useH2opacity):
					#thick case (note that 1.25d-10 = 1/8e9)
					row = row.replace("#KROME_H2opacity", "* min(1.d0, (1.25d-10 * sum(n(1:nmols)))**(-.45))") 
				else:
					#thin case
					row = row.replace("#KROME_H2opacity", "") 

				#replace pragma for dH_cooling
				if(self.useCoolingdH):
					row = row.replace("#KROME_vars","real*8::"+(",".join(dH_varsa))+"\n")
					row = row.replace("#KROME_rates",dH_coe)
					row = row.replace("#KROME_dH_cooling",dH_cool)
				#replace pragma of bremsstrahlung ions
				row = row.replace("#KROME_brem_ions",bms_ions)

				if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"


	##################################################################
	def makeHeating(self):

		reacts = self.reacts
		buildFolder = self.buildFolder
		#*********HEATING****************
		#write header in krome_heating.f90
		print "- writing krome_heating.f90...",

		fh = open("src/krome_heating.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_heating.f90","w")

		#create coefficients, flux, and cooling functions for dH cooling
		i = 0
		dH_varsa = []
		dH_coe = dH_heat = ""
		if(self.useHeatingdH):
			for rea in reacts:
				if(rea.dH!=None and rea.dH>0.e0):
					i += 1 #count heating reactions
					kvar = "k"+str(i) #local variable for coefficient
					dH_varsa.append(kvar) #variables array
					dH_coe += kvar+" = 0.d0\n"
					dH_coe += "if(Tgas."+self.TlimitOpLow+"."+rea.Tmin+" .and. Tgas."+self.TlimitOpHigh+"."+rea.Tmax+") then\n"
					dH_coe += kvar+" = "+rea.krate+"\n" #evaluate reation coefficient
					dH_coe += "end if\n\n" #evaluate reation coefficient
					dH_heat += "heat = heat + "+kvar+"*"+("*".join(["n("+x.fidx+")" for x in rea.reactants])) #evautate heating
					dH_heat += "*"+(str(abs(rea.dH))).replace("e","d")+"\n"


		#build H2 heating according to the rates
		HChem = HChemDust = ""
		sclist = [] 
		if(self.useHeatingChem or self.useCoolingChem):
			RPK = []
			if(self.useHeatingChem):
				RPK.append([["H","H","H"], ["H2","H"], "4.48d0*h2heatfac","H"])
				RPK.append([["H2","H","H"], ["H2","H2"], "4.48d0*h2heatfac","H"])
				RPK.append([["H-","H"], ["H2","E"], "3.53d0*h2heatfac","H"])
				RPK.append([["H2+","H"], ["H2","H+"], "1.83d0*h2heatfac","H"])
			if(self.useCoolingChem):
				RPK.append([["H","E"], ["H+","E","E"], "-13.6d0","C"])
				RPK.append([["HE","E"], ["HE+","E","E"], "-24.6d0","C"])
				RPK.append([["HE+","E"], ["HE++","E","E"], "-79.d0","C"])
				RPK.append([["H2","H"], ["H","H","H"], "-4.48d0","C"])
				RPK.append([["H2","E"], ["H","H","E"], "-4.48d0","C"])
				RPK.append([["H2","H2"], ["H2","H","H"], "-4.48d0","C"])

			Rref = []
			Pref = []
			kref = []
			href = []
			for rpk in RPK:
				Rref.append(sorted(rpk[0]))
				Pref.append(sorted(rpk[1]))
				kref.append(rpk[2])
				href.append(rpk[3])

			for rea in reacts:
				R = sorted([x.name for x in rea.reactants])
				P = sorted([x.name for x in rea.products])
				rmult = ("*".join(["n("+x.fidx+")" for x in rea.reactants]))
				for i in range(len(Rref)):
					if(Rref[i]==R and Pref[i]==P):
						headchem = "!"+rea.verbatim + " ("+("heating" if href[i]=="H"  else "cooling") + ")\n"
						tklim = headchem + "if(Tgas." + rea.TminOp + "."  + rea.Tmin
						tklim += " .and. Tgas." + rea.TmaxOp + "." + rea.Tmax + ") then\n"
						HChem += tklim + "HChem = HChem + k("+str(rea.idx)+") * ("+kref[i] + "*"+rmult+")\n end if\n\n"
						break
			if(self.useDustH2):
				HChemDust += "HChem = HChem + nH2dust * (0.2d0/h2heatfac + 4.2d0)\n"


		#build heating terms for photoionization
		pheatvars = []
		if(self.usePhIoniz):
			for react in reacts:
				phstuff = get_ph_stuff(react)
				if(phstuff==None): continue
				reaname = phstuff["reaname"]
				reag = react.reactants
				fake_opacity = ""
				if(self.useFakeOpacity): fake_opacity = " * exp(-n(" + reag[0].fidx + ") / n0)"
				pheatvars.append("krome_pheat_"+reaname + " * n(" + reag[0].fidx + ")" + fake_opacity)

		#replace pragma with strings built above
		skip = False
		for row in fh:
			if(row.strip() == "#KROME_header"):

				fout.write(get_licence_header(self.version, self.codename))
			else:
				if(row.strip() == "#IFKROME_useHeatingdH" and (not(self.useHeatingdH) or len(dH_varsa)==0)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingCompress" and not(self.useHeatingCompress)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingPhoto" and not(self.useHeatingPhoto)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingChem" and not(self.useHeatingChem) and not(self.useCoolingChem)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(skip): continue

				row = row.replace("#KROME_photo_heating", "photo_heating = " + (" &\n+ ".join(pheatvars)))
				row = row.replace("#KROME_HChem_terms", HChem) #replace chemical heating terms
				row = row.replace("#KROME_HChem_dust", HChemDust) #replace chemical heating for dust
		
				#replace shortcuts for temperature
				if(row.strip() == "#KROME_Tshortcuts"):
					ssc = ""
					for shortcut in sclist:
						ssc += shortcut + "\n"
					row = ssc
				#replace pragma for dH_heating
				if(self.useHeatingdH):
					row = row.replace("#KROME_vars","real*8::"+(",".join(dH_varsa))+"\n")
					row = row.replace("#KROME_rates",dH_coe)
					row = row.replace("#KROME_dH_heating",dH_heat)
			
				if(len(row)==0): continue
				if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()
		print "done!"



	########################################
	def makeODE(self):
		buildFolder = self.buildFolder
		dustArraySize = self.dustArraySize
		dustTypes = self.dustTypes
		nmols = self.nmols
		specs = self.specs
		dnw = self.dnw
		neq = len(specs)
		solver_MF = self.solver_MF
		Tgas_species = self.Tgas_species
		coevarsODE = self.coevarsODE

		#*********ODE****************
		#write parameters in krome_ode.f90
		print "- writing krome_ode.f90...",

		fh = open("src/krome_ode.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_ode.f90","w")

		#string for the function computing dust H2 formation
		dustH2 = "\n"
		if(self.useDustH2):
			iType = 0
			for dType in self.dustTypes:
				ilow = nmols + iType * self.dustArraySize + 1
				iup = nmols + (iType + 1) * self.dustArraySize
				dustH2 += "nH2dust = nH2dust + krome_H2_dust(n(" + str(ilow) + ":" + str(iup) + ")," 
				dustH2 += " krome_dust_T(" + str(ilow-nmols) + ":" + str(iup-nmols) + "), n(:), H2_eps_"+dType+", vgas)\n"
				iType += 1

		#replace pragma with built strings 
		skip = False
		for row in fh:
			srow = row.strip()
			if(srow == "#IFKROME_use_thermo" and (not(self.use_thermo) or not(self.useODEthermo))): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_report" and not(self.doReport)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useDust" and not(self.useDust)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			if(srow == "#KROME_ODE"):
				if(self.use_implicit_RHS):
					fout.write(get_implicit_ode(self.maxnreag, self.maxnprod)+"\n")
				else:
					#add dust ODE and partner specie RHS terms
					if(self.useDust):
						ndust = self.dustArraySize*self.dustTypesSize #number of dust ODEs
						#nmols = len(specs)-4-ndust #number of mols ODEs

						#print first dust to sum gas-dust terms up (e.g. dustC<->gasC)
						for x in dnw[nmols:nmols+ndust]:
							fout.write("\t" + x + "\n")
						fout.write("\n") #print a blank line

						#print sums for different partners
						iType = 0
						for dType in dustTypes:
							ilow = nmols + iType * dustArraySize + 1 #lower index for dust in specs array
							iup = nmols + (iType + 1) * dustArraySize #upper index for dust in specs array
							kpart = "krome_dust_partner_ratio(" + str(ilow-nmols) + ":" + str(iup-nmols) + ")"
							fout.write("\tdSumDust"+dType + " = sum(dn(" + str(ilow) + ":" + str(iup) + ")*"+kpart+")\n")
							iType += 1
						fout.write("\n") #print a blank line

						#print mols ODEs and look for dust partners to add summations
						# and H2 formation on dust and H depletion
						idnw = 0
						for x in dnw[:nmols]:
							for dType in dustTypes:
								if(dType==specs[idnw].name): x += " - dSumDust"+dType
							if(self.useDustH2):
								if("H"==specs[idnw].name): x += " - 2d0*nH2dust"
								if("H2"==specs[idnw].name): x += " + nH2dust"
							fout.write("\t" + x + "\n")
							idnw += 1
						#print other species (CR, PHOTONS, Tgas, dummy)
						for x in dnw[nmols+ndust:]:
							fout.write("\t" + x + "\n")

					#simply write ODEs without dust
					else:
						inw = 0
						for x in dnw:
							#add custom ODE if needed
							if(len(self.customODEs)>0):
								for ode in self.customODEs:
									#build custom ODE
									if(ode[0]==specs[inw].name):
										x = "dn("+str(inw+1)+") = "+ode[1]
										break
							fout.write("\n")
							fout.write("!"+specs[inw].name+"\n")
							fout.write("\t" + x + "\n")
							inw += 1
			elif(srow == "#KROME_calc_Tdust" and self.useDustT):
				getTdust = "nd = " + str(self.dustArraySize) + "\n"
				itype = 0 #dust type index
				#loop in dust types
				for dType in self.dustTypes:
					itype += 1 #increase index
					getTdust += "krome_dust_T(nd*("+str(itype-1)+")+1:nd*"+str(itype) + ") = getTdust("
					getTdust += "dust_opt_Tbb_"+dType+", dust_opt_Em_"+dType+",&\n dust_opt_asize_"+dType+","
					getTdust += " dust_opt_nu_"+dType+", dust_opt_Qabs_"+dType+", n(:))\n"
					getTdust = getTdust.replace("nd*(0)+1","1").replace("nd*1","nd").replace("nd*(1)","nd")
				fout.write(getTdust+"\n")

			elif(srow == "#KROME_initcoevars"):
				if(len(coevarsODE)==0): continue
				kvars = "real*8::"+(",".join([x for x in coevarsODE.keys()]))
				fout.write(kvars+"\n")
			elif(srow == "#KROME_coevars"):
				if(len(coevarsODE)==0): continue
				klist = [[k+" = "+v[1]+"\n",v[0]] for k,v in coevarsODE.iteritems()] #this mess is to sort dict
				klist = sorted(klist, key=lambda x: x[1])
				fout.write("".join([x[0] for x in klist]))

			elif(srow == "#KROME_dustSumVariables" and self.useDust):
				#add partner sum dust variable declarations
				dustSumVar = []
				for dType in dustTypes:
					dustSumVar.append("dSumDust"+dType)
				fout.write("\t real*8::" + (",".join(dustSumVar)) + "\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			elif(srow == "#KROME_implicit_variables"):
				fout.write("real*8::rr\n")
				ris = (",".join(["r"+str(i+1) for i in range(self.maxnreag)]))
				pis = (",".join(["p"+str(i+1) for i in range(self.maxnprod)]))
				rpis = ",".join([x for x in ["i",ris,pis] if(len(x)>0)])
				fout.write("integer::"+rpis+"\n")
			elif(srow == "#KROME_report_flux"):
				report_flux = ("*".join(["n(arr_r"+str(j+1)+"(i))" for j in range(self.maxnreag)]))
		 		fout.write("write(fnum,'(I5,E12.3e3,a2,a50)') i,k(i)*"+report_flux+",'',rnames(i)\n")

			elif(srow == "#KROME_odeConstant" and self.useODEConstant):
				fout.write("dn(:) = dn(:) "+self.ODEConstant+" \n") #add the string contains an ODE expression
			elif(srow == "#KROME_odeDust"):
				fout.write("dn(:) =  krome_dust \n")
			elif(srow == "#KROME_dust_H2"):
				fout.write(dustH2+"\n")
			elif(srow == "#KROME_JAC_PD"):
				if(solver_MF==222):
					fout.write("\n")
				else:
					#build the Jacobian as J(i,j) = df_i/dx_j
					for i in range(neq):
						if(i==0): fout.write("if(j=="+str(i+1)+") then\n")
						if(i>0): fout.write("elseif(j=="+str(i+1)+") then\n")
						if(i!=Tgas_species.idx-1):
							spdj = ""
							has_pdj = False
							for j in range(neq):
								if(self.jsparse[j][i]==1):
									has_pdj = True
									spdj += ("\t" + self.jac[j][i] + "\n")
							#if(has_pdj): fout.write("k(:) = coe_tab(n(:))\n")
							fout.write(spdj)
						else:
							jacT = "!use fex to compute temperature-dependent Jacobian\n"
							if(self.deltajacMode=="RELATIVE"):
								jacT += "dnn = n(idx_Tgas)*"+str(self.deltajac)+"\n"
							elif(self.deltajacMode=="ABSOLUTE"):
								jacT += "dnn = "+str(self.deltajac)+"\n"
							else:
								die("ERROR: unknown deltajacMode! "+self.deltajacMode)
							jacT += """nn(:) = n(:)
								nn(idx_Tgas) = n(idx_Tgas) + dnn
								call fex(neq,tt,nn(:),dn(:))
								do i=1,neq-1
								  pdj(i) = dn(i) / dnn
								end do"""
							if(not(self.use_thermo)): jacT = ""
							fout.write("\t" + jacT.replace("\t","") + "\n")
					fout.write("end if\n")
			elif(srow == "#KROME_gamma"):
				#computes the adiabatic index if needed or uses a user-defined expression 
				if(self.typeGamma=="DEFAULT"):
					gamma = "1.66666666667d0" #default gamma is 5/3 (atomic)
				elif(self.typeGamma=="FULL"):
					#build gamma following (Grassi+2010,MerlinPhDTheis)
					specsGamma = ["H","HE","E","H2"] #species used in the gamma
					gammaNm = [] #numerator monoatomics
					gammaDm = [] #denominator monoatomics
					gammaNb = [] #numerator diatomics
					gammaDb = [] #denominator diatomics
					#loop on species
					for mol in specs:
						if(mol.name.upper() in specsGamma):
							nfidx = "n("+mol.fidx+")" #number density
							if(mol.is_atom): #monoatomic
								gammaNm.append(nfidx)
								gammaDm.append(nfidx)
							else: #diatomic
								gammaNb.append(nfidx)
								gammaDb.append(nfidx)
					#build fortran expression for gamma
					gammaN = "(5.d0*("+(" + ".join(gammaNm)) + ") + 7.d0*("+(" + ".join(gammaNb)) + "))"
					gammaD = "(3.d0*("+(" + ".join(gammaDm)) + ") + 5.d0*("+(" + ".join(gammaDb)) + "))"
					gamma = gammaN + " / &\n" +gammaD
				else: 
					#user-defined gamma
					gamma = self.typeGamma

				#write gamma
				fout.write("krome_gamma = " + gamma + "\n")
		
			else:
				srow = row.strip()
				if(len(srow)>0):
					if(srow[0]!="#"): fout.write(row)
				else:
					fout.write(row)
		if(not(self.buildCompact)):
			fout.close()
		print "done!"
	################################
	def makeUser(self):

		reacts = self.reacts
		nmols = self.nmols
		dustArraySize = self.dustArraySize
		dustTypesSize = self.dustTypesSize
		buildFolder = self.buildFolder
		specs = self.specs
		#*********USER****************
		#write parameters in krome_user.f90
		print "- writing krome_user.f90...",
		fh = open("src/krome_user.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_user.f90","w")

		skip = False

		#solar abundances from Anders+Grevesse 1989
		solar = {
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

		scaleZ = []
		#looks for H to rescale the metallicity otherwise skips
		has_H = False
		for mols in specs:
			if(mols.name=="H"): 
				has_H = True
				break

		#creates the metallicity rescaling subroutine
		for (k,v) in solar.iteritems():
			if(not(has_H)): break #skip routine if H is not present
			for mols in specs:
				if(k.upper()==mols.name.upper()):
					scaleZ.append("n("+mols.fidx+") = max(n(idx_H) * 1d1**(Z+log10("+str(v)+")), 1d-40)")

		for row in fh:

			srow = row.strip()
			if(srow == "#IFKROME_use_cooling" and not(self.use_cooling)): skip = True
			if(srow == "#ENDIFKROME"): skip = False
			if(skip): continue

			if(srow == "#KROME_species"):
				for x in specs:
					fout.write("\tinteger,parameter::" + "KROME_"+x.fidx + " = " + str(x.idx) +"\t!"+x.name+"\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			elif(srow == "#KROME_zero_electrons"):
				#check if electron exists
				for x in specs:
					if(x.name=="E"):
						fout.write("x(idx_e) = 0.d0\n")
						break
			elif(srow == "#KROME_constant_list"):
				const = ""
				constants = self.constantList
				newc = []
				for i in range(len(constants)):
					x = constants[i]
					for j in range(i):
						y = constants[j]
						x[1] = x[1].replace(y[0],"krome_"+y[0])
					newc.append(x)

				for x in newc:
					const += "real*8,parameter::krome_" + x[0] + " = " + x[1] + " !" + x[2] + "\n"
				fout.write(const)
			elif(srow == "#KROME_common_alias"):
				fout.write("\tinteger,parameter::krome_nrea=" + str(len(reacts)) + "\n")
				fout.write("\tinteger,parameter::krome_nmols=" + str(nmols) + "\n")
				fout.write("\tinteger,parameter::krome_nspec=" + str(len(specs)) + "\n")
				fout.write("\tinteger,parameter::krome_ndust=" + str(dustArraySize*dustTypesSize) + "\n")
				fout.write("\tinteger,parameter::krome_ndustTypes=" + str(dustTypesSize) + "\n")
			elif(srow == "#KROME_scaleZ"):
				fout.write(("\n".join(scaleZ))+"\n")
			else:
				if(len(srow)>0):
					if(srow[0]!="#"): fout.write(row)
				else:
					fout.write(row)
		if(not(self.buildCompact)):
			fout.close()
		print "done!"

	####################################
	def makeReduction(self):
		buildFolder = self.buildFolder
		#********* REDUCTION ****************
		#WARNING: this part is not supported and its use is discouraged
		print "- writing krome_reduction.f90...",
		fh = open("src/krome_reduction.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_reduction.f90","w")

		skip = False
		for row in fh:
			srow = row.strip()
			if(srow == "#IFKROME_useTopology" and not(self.useTopology)): skip = True
			if(srow == "#ENDIFKROME"): skip = False
			if(srow == "#IFKROME_useFlux" and not(self.useFlux)): skip = True
			if(srow == "#ENDIFKROME"): skip = False
			if(srow == "#IFKROME_useReduction" and not(self.useTopology) and not(self.useFlux)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue
			if(srow == "#KROME_rvars"):
				if(self.maxnreag>0):
					fout.write("integer::"+(",".join(["r"+str(j+1) for j in range(self.maxnreag)]))+"\n")
			if(srow == "#KROME_arrs"):
				for j in range(self.maxnreag):
					fout.write("r"+str(j+1)+" = arr_r"+str(j+1)+"(i)\n")
			if(srow == "#KROME_arr_flux"):
				print self.maxnreag
				if(self.maxnreag>0):
					fout.write("arr_flux(i) = k(i)*"+("*".join(["n(r"+str(j+1)+")" for j in range(self.maxnreag)]))+"\n")

			if(row[0]!="#"): fout.write(row)
		if(not(self.buildCompact)):
			fout.close()

		print "done!"

	###############################
	def makeMain(self):
		buildFolder = self.buildFolder
		dustArraySize = self.dustArraySize
		dustTypes = self.dustTypes
		#*********MAIN****************
		#write WORKS arrays and IAC/JAC in krome.f90
		print "- writing krome.f90...",
		if(self.useDvodeF90):
			fh = open("src/kromeF90.f90")
		else:
			fh = open("src/krome.f90")

		ATOL = self.ATOL
		RTOL = self.RTOL
		if(is_number(ATOL)): ATOL = '%e' % ATOL
		if(is_number(RTOL)): RTOL = '%e' % RTOL
		ATOL = ATOL.replace("e","d")
		RTOL = RTOL.replace("e","d")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome.f90","w")

		skip = False
		for row in fh:
			srow = row.strip()
			if(srow == "#IFKROME_useX" and not(self.useX)): skip = True
			if(srow == "#ELSEKROME" and not(self.useX)): skip = False
			if(srow == "#ELSEKROME" and self.useX): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useTabs" and not(self.useTabs)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useFlux" and not(self.useFlux)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_report" and not(self.doReport)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useTopology" and not(self.useTopology)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_check_mass_conservation" and not(self.checkConserv)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useDust" and not(self.useDust)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useEquilibrium" and not(self.useEquilibrium)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(self.useDust):
				row = row.replace("#KROME_dust_arguments",",xdust")
			else:
				row = row.replace("#KROME_dust_arguments","")

			row = row.replace("#KROME_ATOL",str(ATOL))
			row = row.replace("#KROME_RTOL",str(RTOL))

			if(skip): continue

			if(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename))
			elif(srow == "#KROME_custom_ATOL"):
				#add custom atols
				if(len(self.atols)>0):
					for x in self.atols:
						#check for species in species list
						for y in self.specs:
							#one can use either the name or the idx name (e.g. H+ or idx_Hp)
							if(x[0]==y.name or x[0]==y.fidx):
								fout.write("atol("+y.fidx+") = "+format_double(x[1])+"\n")
								break
			elif(srow == "#KROME_custom_RTOL"):
				#add custom rtols
				if(len(self.rtols)>0):
					for x in self.rtols:
						#check for species in species list
						for y in self.specs:
							#one can use either the name or the idx name (e.g. H+ or idx_Hp)
							if(x[0]==y.name or x[0]==y.fidx):
								fout.write("rtol("+y.fidx+") = "+format_double(x[1])+"\n")
								break
			elif(srow == "#KROME_rwork_array"):
				fout.write("\treal*8::rwork("+str(self.lrw)+")\n")
			elif(srow == "#KROME_iwork_array"):
				fout.write("\tinteger::iwork("+str(30+len(self.ia)+len(self.ja))+")\n")
			elif(srow == "#KROME_init_IAC"):
				fout.write("\t"+self.iaf+"\n")
			elif(srow == "#KROME_init_JAC"):
				fout.write("\t"+self.jaf+"\n")
			elif(srow == "#KROME_maxord" and self.maxord!=0):
				fout.write("iopt = 1 !activate optional inputs\n")
				fout.write("IWORK(5) = "+str(self.maxord)+" !maximum integration order\n")
			elif(srow == "#KROME_MF"):
				if(self.solver_MF==21):
					fout.write("\tMF=21\n")
				elif(self.solver_MF==222):
					fout.write("\tMF=222\n")
			else:
				if(row[0]!="#"): fout.write(row)
		if(not(self.buildCompact)):
			fout.close()
		print "done!"

	###################################
	def makeReport(self):
		specs = self.specs
		neq = len(specs)
		#*********REPORT.gps****************
		#write gnuplot script to plot abundances evolution
		if(self.doReport):
			fout = open(self.buildFolder+"report.gps","w")
			fout.write("#plot KROME report\n")
			fout.write("reset\n")
			fout.write("set logscale\n")
			fout.write("set autoscale\n")
			fout.write("set xlabel \"log(time/s)\"\n")
			fout.write("plot 'fort.98' u 1:2 w l t \""+specs[0].name+"\",\\\n")
			for i in range(neq-2):
				fout.write("'' u 1:"+str(i+3)+" w l t \""+specs[i+1].name+"\",\\\n")
			fout.write("'' u 1:"+str(neq+2)+" w l t \""+specs[neq-1].name+"\"\n")
			fout.close()

	############################################
	def copyOthers(self):
		buildFolder = self.buildFolder
		test_name = self.test_name
		#copy other files to build
		print "- copying others...",

		print " copying optical data for dust..."
		if(self.useDust):
			shutil.copyfile("data/optC.dat", buildFolder+"optC.dat")
			shutil.copyfile("data/optSi.dat", buildFolder+"optSi.dat")

		#copy static files to build
		if(self.is_test):
			print "- copying test to /build..."
			mypath = "tests/"+test_name
			files = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
			for fdir in files:
				shutil.copyfile("tests/"+test_name+"/"+fdir, buildFolder+"/"+fdir)
				print "- copying "+fdir+" to "+buildFolder
			if(self.useDvodeF90):
				shutil.copyfile("tests/MakefileF90", buildFolder+"Makefile")
			elif(self.buildCompact):
				shutil.copyfile("tests/"+test_name+"/MakefileCompact", buildFolder+"Makefile")
			else:
				if(self.pedanticMakefile):
					shutil.copyfile("tests/Makefile_pedantic", buildFolder+"Makefile")
				else:
					shutil.copyfile("tests/"+test_name+"/Makefile", buildFolder+"Makefile")

			#test_file = "tests/"+test_name+"/test.f90"
			#plot_file = "tests/"+test_name+"/plot.gps"

			#chech if test file exists
			#if(not(os.path.isfile(filename))): die("ERROR: Test file \""+test_file+"\" doesn't exist!")

			#shutil.copyfile(test_file, buildFolder+"test.f90")
			#if(self.has_plot): shutil.copyfile(plot_file, buildFolder+"plot.gps")
			#print " done!"

		#copy solver files to build
		print "- copying solver to /build...",
		if(self.useDvodeF90):
			shutil.copyfile("solver/dvode_f90_m.f90", buildFolder+"dvode_f90_m.f90")
		else:
			shutil.copyfile("solver/opkdmain.f", buildFolder+"opkdmain.f")
			shutil.copyfile("solver/opkda1.f", buildFolder+"opkda1.f")
			shutil.copyfile("solver/opkda2.f", buildFolder+"opkda2.f")
		print " done!"


	#######################################################
	def indent(self):
		buildFolder = self.buildFolder
		print "Indenting...",
		if(self.doIndent):
			if(self.buildCompact):
				indentF90(buildFolder+"krome_all.f90")
			else:
				indentF90(buildFolder+"krome_commons.f90")
				indentF90(buildFolder+"krome_constants.f90")
				indentF90(buildFolder+"krome_cooling.f90")
				indentF90(buildFolder+"krome_heating.f90")
				indentF90(buildFolder+"krome_dust.f90")
				indentF90(buildFolder+"krome.f90")
				indentF90(buildFolder+"krome_ode.f90")
				indentF90(buildFolder+"krome_photo.f90")
				indentF90(buildFolder+"krome_reduction.f90")
				indentF90(buildFolder+"krome_subs.f90")
				indentF90(buildFolder+"krome_tabs.f90")
				indentF90(buildFolder+"krome_user.f90")

		print "done!"
	#########################################
	#copy fsrc to fout file and replace the list in pragmas with the list in repls
	def replacein(self,fsrc,fout,pragmas,repls,trim=True):
		fh = open(fsrc,"rb")
		fw = open(fout,"w")
		if(len(pragmas)!=len(repls)):
			print "ERROR: in replacein len(pragmas)!=len(repls)"
			sys.exit()
		for row in fh:
			srow = row
			if(trim): srow = row.strip()
			#replace only with non-empty lists
			if(len(pragmas)>0):
				for i in range(len(pragmas)):
					x = pragmas[i]
					y = str(repls[i])
					srow = srow.replace(x,y)
			if(trim): 
				fw.write(srow+"\n")
			else:
				fw.write(srow)
		fh.close()
		fw.close()


				

	#########################################
	def ramses_patch(self):
		pfold = "patches/ramses/"
		ramsesFolder = self.buildFolder+"krome_ramses_patch/" 
		buildFolder = self.buildFolder 
		if not os.path.exists(ramsesFolder): os.makedirs(ramsesFolder)
		specs = self.specs

		#some initial abundances. if not found set default
		ndef = {"H": 7.5615e-1,
			"E": 4.4983e-8,
			"H+": 8.1967e-5,
			"HE": 2.4375e-1,
			"H2": 1.5123e-6,
			"default":1e-20
		}

		excl = ["CR","g","Tgas","dummy"] #avoid specials

		#count species excluding excl
		chemCount = 0
		for x in specs:
			if(x.name in excl): continue
			chemCount += 1

		#amr_parameters
		fname = "amr_parameters.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,["aaa"],["aaa"])
		indentF90(ramsesFolder+fname)

		#condinit
		cheminit = " q(1:nn,ndim+3) = 200.d0     !Set temperature in K\n"
		ichem = 0
		fname = "condinit.f90"
		for x in specs:
			if(x.name in excl): continue
			ichem += 1
			#check if species is ini init array
			if(x.name in ndef):
				sdef = str(ndef[x.name]) #default value from array
			else:
				sdef = str(ndef["default"]) #default values if not present in array
			cheminit += "q(1:nn,ndim+3+"+str(ichem)+")  = "+sdef+"  !"+x.name+"\n"
		self.replacein(pfold+fname,ramsesFolder+fname,["#KROME_init_chem"],[cheminit])
		indentF90(ramsesFolder+fname)

		#cooling_fine
		updateueq = scaleueq = bkscaleueq = bkupdateueq = ""
		ichem = 0
		fname = "cooling_fine.f90"
		for x in specs:
			ichem += 1
			if(not(x.name in excl)):
				updateueq += "unoneq("+str(ichem)+") = uold(ind_leaf(i),ndim+3+"+str(ichem)+") !"+x.name+"\n"
				if(x.mass>0e0): scaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*scale_d/"+str(x.mass)+" !"+x.name+"\n"
				bkscaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*"+str(x.mass)+"/scale_d !"+x.name+"\n"
				bkupdateueq += "uold(ind_leaf(i),ndim+3+"+str(ichem)+") = unoneq("+str(ichem)+")\n"
		org = ["#KROME_update_unoneq","#KROME_scale_unoneq","#KROME_backscale_unoneq","#KROME_backupdate_unoneq"]
		new = [updateueq, scaleueq, bkscaleueq, bkupdateueq]
		self.replacein(pfold+fname,ramsesFolder+fname,org,new)
		indentF90(ramsesFolder+fname)

		#cooling_module
		fname = "cooling_module.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		indentF90(ramsesFolder+fname)

		#hydro_parameters
		fname = "hydro_parameters.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		#indentF90(ramsesFolder+fname)

		#init_flow_fine
		fname = "init_flow_fine.f90"
		init_array = "if(ivar==ndim+3)  init_array = 1.356d-2/aexp**2 ! T in K\n"
		ichem = 0
		for x in specs:
			if(x.name in excl): continue
			ichem += 1
			#check if species
			if(x.name in ndef):
				sdef = str(ndef[x.name]) #default value from array
			else:
				sdef = "1d-20" #default values if not present in array
			init_array += "if(ivar==ndim+3+"+str(ichem)+")  init_array = "+sdef+"  !"+x.name+"\n"
		self.replacein(pfold+fname,ramsesFolder+fname,["#KROME_init_array"],[init_array])
		indentF90(ramsesFolder+fname)
		
		#output_hydro
		fname = "output_hydro.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		indentF90(ramsesFolder+fname)

		#read_hydro_params
		fname = "read_hydro_params.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		indentF90(ramsesFolder+fname)

		#Makefile
		fname = "Makefile"
		#note that makefile will be copied in the build folder
		self.replacein(pfold+fname,buildFolder+fname,["#KROME_nvar"],["#this must be NDIM+"+str(chemCount+1)], False)

		#move the krome files into the ramses patch folder
		shutil.move(buildFolder+"krome_all.f90", ramsesFolder+"krome_all.f90")
		shutil.move(buildFolder+"krome_user_commons.f90", ramsesFolder+"krome_user_commons.f90")
		shutil.move(buildFolder+"opkda1.f", ramsesFolder+"opkda1.f")
		shutil.move(buildFolder+"opkda2.f", ramsesFolder+"opkda2.f")
		shutil.move(buildFolder+"opkdmain.f", ramsesFolder+"opkdmain.f")



	###########################################
	def flash_patch(self):
		return
	
	###########################################
	def enzo_patch(self):
		return

	############################################
	def patches(self):
		if(self.doFlash): self.flash_patch()
		if(self.doRamses): self.ramses_patch()
		if(self.doEnzo): self.enzo_patch()
		return
		if(self.doFlash or self.doRamses or self.doEnzo):		
			print "**********************************************************"
			print "  We're sorry, but in this version of KROME"
			print "  the patch builder is available only on request."
			print "  If you are interested please contact the krome"
			print "  developers to obtain the beta versions of the patches."
			print "  We will give all the support you need!"
			print
			print "  email       : krome@kromepackage.org"
			print "  www         : http://kromepackage.org/"
			print "  mailing list: https://groups.google.com/forum/#!forum/kromeusers"
			print "  bitbucket   : https://bitbucket.org/krome/krome_stable"
			print "**********************************************************"
			sys.exit()

	#########################################
	def final_report(self):

		buildFolder = self.buildFolder
		reacts = self.reacts
		useX = self.useX
		nmols = self.nmols
		print
		print "You'll find the necessary files in "+buildFolder
		print "Example call to the solver in "+buildFolder+"test.f90"
		print "Example Makefile in "+buildFolder+"Makefile"

		#check for large reaction set
		if(len(reacts)>500 and not(self.use_implicit_RHS)):
			print
			print "WARNING: "+str(len(reacts))+" reactions found! Using implicit RHS (option -iRHS) could be more efficient"
			print "  and also allows faster compilation."
			a = raw_input("Any key to continue q to quit... ")
			if(a=="q"): print sys.exit()

		print
		#IF NOT TEST
		if(not(self.is_test)):
			#PATCHES DO NOT NEED MAKEFILE AND TEST.F90
			if(not(self.doFlash or self.doRamses or self.doEnzo)):
				#TODO: add description in case of dust
				print "Call KROME in your code as:"
				if(useX):
					print "    call krome(x(:), gas_density, gas_temperature, time_step)"
				else:
					print "    call krome(x(:), gas_temperature, time_step)"
				print "where:" 
				print " x(:) is a real*8 array of size "+str(nmols)+(" of the mass fractions" if useX else " of number densities [1/cm3]")
				if(useX): print " gas_density  is the gas density in [g/cm3]"
				print " gas_temperature is the gas temperature in [K]"
				print " time_step is the integration time-step in [s]"
			
				fout = open(buildFolder+"test.f90","w")
				fout.write(get_example(nmols,useX))
				fout.close()
				indentF90(buildFolder+"test.f90")
				if(self.buildCompact):
					shutil.copyfile("tests/MakefileCompact", buildFolder+"Makefile")
				else:
					if(self.pedanticMakefile):
						shutil.copyfile("tests/Makefile_pedantic", buildFolder+"Makefile")
					else:
						shutil.copyfile("tests/Makefile", buildFolder+"Makefile")

		#IF IT IS A TEST
		else:
			print "This is a test. To run it just type:"
			print "> cd build/"
			print "> make"
			print "> ./krome"
		print
		print "Everything done, goodbye!"
