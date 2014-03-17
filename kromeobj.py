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
# University of Rome \"Sapienza\".
#
# Co-developer Stefano Bovino
# sbovino@astro.physik.uni-goettingen.de
# Institut fuer Astrophysik, Goettingen.
#
# Others (alphabetically): F.A. Gianturco, J.Prieto,
# D.R.G. Schleicher, D. Seifried, E. Simoncini , E. Tognelli
#
#
#KROME is provided \"as it is\", without any warranty. 
# The Authors assume no liability for any damages of any kind 
# (direct or indirect damages, contractual or non-contractual 
# damages, pecuniary or non-pecuniary damages), directly or 
# indirectly derived or arising from the correct or incorrect 
# usage of KROME, in any possible environment, or arising from 
# the impossibility to use, fully or partially, the software, 
# or any bug or malefunction.
# Such exclusion of liability expressly includes any damages 
# including the loss of data of any kind (including personal data)

# THIS FILE CONTAINS THE KROME CLASS

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
	useCoolingCompton = useCoolingExpansion = useH2opacity = useShieldingDB96 = useShieldingWG11 = useCoolingCIE = useCoolingDISS  = False
	#useCoolingZC = useCoolingZCp = useCoolingZSi = useCoolingZSip = useCoolingZO = useCoolingZOp = useCoolingZFe = useCoolingZFep = False
	useReverse = useCustomCoe = useODEConstant = cleanBuild = usePlainIsotopes = useDust = use_thermo = useStars = useNuclearMult = False
	usePhIoniz = useHeatingCompress = useHeatingPhoto = useHeatingChem = useDecoupled = useCoolingdH = useHeatingdH = useCoolingChem = False
	pedanticMakefile = useFakeOpacity = useConserve = useConserveE = noExample = False
	useX = has_plot = doIndent = useTlimits = useODEthermo = safe = doJacobian = True
	useDustGrowth = useDustSputter = useDustH2 = useDustT = checkThermochem = needLAPACK = False
	doRamses = doRamsesTH = doFlash = doEnzo = wrapC = mergeTlimits = shortHead = isdry = useIERR = checkReverse = False
	humanFlux = True
	typeGamma = "DEFAULT"
	test_name = "default"
	is_test = False
	TlimitOpLow = "GE"
	TlimitOpHigh = "LT"
	customCoeFunction = "[CUSTOM COE FUNCTION NOT SET!]"
	buildFolder = "build/"
	srcFolder = "src/"
	TminAuto = 1e99
	TmaxAuto = 0e0
	checkMode = "ALL" #conservation check mode (ALL | [CHARGE],[MASS]| NONE)
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
	jaca = [] #unrolled sparse jacobian
	customODEs = [] #custom ODEs
	nrea = 0 #number of reactions
	full_cool = vars_cool = ""
	coolZ_functions = []
	coolZ_rates = []
	coolZ_vars_cool = []
	coolZ_nkrates = 0
	zcoolants = [] #list of cooling read from file (flag name, e.g CII)
	Zcools = [] #list of cooling read from file (species name, e.g. C+)
	ramses_offset = 3 #offset in the array for ramses
	coolFile = "data/coolZ.dat"
	version = "14.03"
	codename = "Beastie Boyle"

	#########################################
	def checkPrereq(self):
		#check for argparse module
		try:
			import argparse
		except:
			print "ERROR: you need installed argparse!"
			print "You can obtain it by typing (ubuntu users):"
			print " apt-get install python-setuptools"
			print " easy_install argparse"
			print ""
			print "more details here:"
			print " https://pypi.python.org/pypi/argparse"
			sys.exit()

		#check python version
		ver = sys.version_info
		aver = list(ver)
		sver = (".".join([str(x) for x in aver[:3]]))
		if(not(aver[0]>=2 and aver[1]>=5)):
			print "ERROR: your version of Python ("+sver+") is not supported by KROME!"
			print " KROME needs at least Python 2.5!"
			sys.exit()

		#check necessary files
		fles = get_file_list() #get the list of necessary files
		for fle in fles:
			if(os.path.isdir(fle)): continue #do not check folders
			if(not(os.path.isfile(fle))):
				print "************************************************"
				print "WARNING: the file "+fle+" is missing!"
				print "Do you want to proceed anyway?"
				print "************************************************"
				a = raw_input("Any key to ignore q to quit... ")
				if(a=="q"): print sys.exit()

	#########################################
	def init_argparser(self):

		tests = ", ".join(os.walk('tests').next()[1])
		self.parser = argparse.ArgumentParser(description="KROME a package for astrochemistry")
		self.parser.add_argument("-ATOL", help="set solver absolute tolerance to the float or double value ATOL, e.g. -atol 1d-40\
			Default is ATOL=1d-20, see also -RTOL and -customATOL")
		self.parser.add_argument("-C", action="store_true", help="create a simple C wrapper")
		self.parser.add_argument("-compact", action="store_true", help="creates a single fortran file with all the modules instead of\
			various file with the different modules. Solver files remain stand-alone (see example make in test/MakefileCompact)")
		self.parser.add_argument("-checkConserv", action="store_true", help="check mass conservation during integration (slower)")
		self.parser.add_argument("-checkReverse", action="store_true", help="check network for reverse reactions. Write warning on\
			screen if any.")
		self.parser.add_argument("-checkThermochem", action="store_true", help="print a warning when thermochemistry data are not found\
			for a given species.")
		self.parser.add_argument("-clean", action="store_true", help="clean all in /build (including krome_user_commons.f90 that\
			is normally kept by default) before creating new f90 files.")
		self.parser.add_argument("-compressFluxes", action="store_true", help="in the ODE fluxes are stored in a single variable")
		self.parser.add_argument("-conserve", action="store_true", help="conserves the species total number and charge global\
			neutrality. Works with some limitations, please read the manual.")
		self.parser.add_argument("-conserveE", action="store_true", help="conserves the charge global neutrality only.")
		self.parser.add_argument("-coolFile", metavar='FILENAME', help="select the filename to be used to load external cooling. See\
			also tools/lamda2.py script for a LAMDA<->KROME converter. Default FILENAME is data/coolZ.dat, which contains\
			fine-strucutre atomic metal cooling for C,O,Si,Fe, and their first ions.")
		self.parser.add_argument("-cooling", metavar='TERMS', help="cooling options, TERMS can be ATOMIC, H2, HD, Z, DH, DUST, H2GP98,\
			COMPTON, EXPANSION, CIE, DISS, CI, CII, SiI, SiII, OI, OII, FeI, FeII, CHEM (e.g. -cooling=ATOMIC,CII,OI,FeI). Note that further\
			cooling options can be added when reading cooling function from file")
		self.parser.add_argument("-customATOL", help="file with the list of the individual ATOLs in the form SPECIES ATOL in each line,\
			e.g. H2 1d-20, see also -ATOL", metavar="filename")
		self.parser.add_argument("-customODE", help="file with the list of custom ODEs", metavar="FILENAME")
		self.parser.add_argument("-customRTOL", help="file with the list of the individual RTOLs in the form SPECIES RTOL in each line,\
			e.g. H3+ 1d-4, see also -RTOL", metavar="filename")
		self.parser.add_argument("-dry", action="store_true", help="dry pre-compilation: does not write anything in the build direactory")
		self.parser.add_argument("-dust", help="include dust ODE using N bins for each TYPE, e.g. -dust 10,C,Si set 10 dust carbon\
			bins and 10 dust silicon dust bins. Note: requires a call to the krome_init_dust subroutine.\
			See -test=dust for an example.")
		self.parser.add_argument("-dustOptions", help="activate dust options: (GROWTH) dust growth, (SPUTTER) sputtering, (H2) molecular\
			hydrogen formation on dust, and (T) dust temperature. The last option provide a template for the FEX routine.",\
			metavar="OPTIONS")
		self.parser.add_argument("-enzo", action="store_true", help="create patches for ENZO")
		self.parser.add_argument("-flash", action="store_true", help="create patches for FLASH")
		self.parser.add_argument("-forceMF21", action="store_true", help="force explicit sparsity and Jacobian")
		self.parser.add_argument("-forceMF222", action="store_true", help="force internal-generated sparsity and Jacobian")
		self.parser.add_argument("-forceRWORK", help="force the size of RWORK to N", metavar="N")
		self.parser.add_argument("-gamma",help="define the adiabatic index according to OPTION that can be FULL for employing Grassi et al.\
			2011, i.e. a density dependent but temperature independent adiabatic index, VIB to keep into account the vibrational\
                        paritition function, ROT to keep into account the rotational partition function, or EXACT to evaluate the\
                        adiabatic index accurately taking into account both contributions. Finally a custom F90 expression e.g. -gamma=\"1d0\"\
			can also be used. Default value is 5/3.",metavar="OPTION")
		self.parser.add_argument("-heating", metavar='TERMS', help="heating options, TERMS can be COMPRESS, PHOTO, CHEM, DH")
		self.parser.add_argument("-ierr", action="store_true", help="same as -useIERR")
		self.parser.add_argument("-iRHS", action="store_true", help="implicit loop-based RHS (suggested for large systems)")
		self.parser.add_argument("-maxord", help="max order of the BDF solver. Default (and maximum values) is 5.")
		self.parser.add_argument("-mergeTlimits", action="store_true", help="use the same reaction index for equivalent\
			reactions (same reactants and products) that have different temperature limits")
		self.parser.add_argument("-n", help="reaction network file", metavar='FILENAME')
		self.parser.add_argument("-nochargeCheck", action="store_true", help="skip reaction charge check")
		self.parser.add_argument("-noCheck", action="store_true", help="skip reaction charge and mass check. Equivalent to\
			-nomassCheck -nochargeCheck options.")
		self.parser.add_argument("-noExample", action="store_true", help="do not write test.f90 and Makefile in the build directory")
		self.parser.add_argument("-nomassCheck", action="store_true", help="skip reaction mass check")
		self.parser.add_argument("-noTlimits", action="store_true", help="ignore rate coefficient temperature limits.")
		self.parser.add_argument("-nuclearMult", action="store_true", help="keep into account reactants multeplicity, and modify\
			fluxes according to this. Intended for nuclear networks.")
		self.parser.add_argument("-options", metavar="filename", help="read the options from a file instead of command line\
			(in principle you can use both)")
		self.parser.add_argument("-pedantic", action="store_true", help="uses a pedantic Makefile (debug purposes)")
		self.parser.add_argument("-project", help="build everything in a folder called build_NAME instead of building all in the\
			default build folder. It also creates a NAME.kpj file with the krome input used.",metavar="NAME")
		self.parser.add_argument("-quote", action="store_true", help="print a citation and exit")
		self.parser.add_argument("-ramses", action="store_true", help="create patches for RAMSES, see also -enzo and -flash")
		self.parser.add_argument("-ramsesOffset", metavar="offset", help="add an offset to the array of the passive scalar. The\
			default is 3.")
		self.parser.add_argument("-ramsesTH", action="store_true", help="create patches for RAMSES_TH. This is a private version\
			and probably does not fix your needs.")
		self.parser.add_argument("-report", action="store_true", help="generate report file in the main call to krome as\
			KROME_ERROR_REPORT and when calling the fex as KROME_ODE_REPORT. It also stores abundances evolution in fex as \
			fort.98, and prepares a report.gps gnuplot script file to plot evolutions callable in gnuplot with load \
			'report.gps'. Warning: it slows the whole system!")
		self.parser.add_argument("-reverse", action="store_true", help="create reverse reaction from the current network\
			using NASA polynomials.")
		self.parser.add_argument("-RTOL", help="set solver relative tolerance to the float double value RTOL, e.g.\
			-RTOL 1e-5 Default is RTOL=1d-4, see also -ATOL and -customRTOL")
		self.parser.add_argument("-sh", action="store_true", help="write a shorter header in the f90 files")
                self.parser.add_argument("-shielding", metavar="TYPE", help="use H2 self-shielding, TYPE can be DB96 for Draine+Bertoldi 1996,\
                        WG11 for the more accurate Wolcott+Greene 2011")
		self.parser.add_argument("-skipDup", action="store_true", help="skip duplicate reactions")
		self.parser.add_argument("-skipJacobian", action="store_true", help="do not write Jacobian in krome_ode.f90 file. Useful\
			to reduce compilation time when Jacobian is not needed (MF=222).")
		self.parser.add_argument("-skipODEthermo", action="store_true", help="do not compute dT/dt in the ODE RHS function (fex)")
		self.parser.add_argument("-source", metavar="folder", help="use FOLDER as source directory")
		self.parser.add_argument("-stars", action="store_true", help="use star module for nuclear reactions. NOTE: krome_stars\
			module required in the Makefile")
		self.parser.add_argument("-test",help=("Create a test model in /build. TEST can be: "+tests+"."))
		self.parser.add_argument("-Tlimit", metavar="opLow,opHigh", help="set the operators for all the reaction temperature limits\
			where opLow is the operator for the first temperature value in the reaction file, and opHigh is for the second one. e.g.\
			if the T limits for a given reaction are 10. and 1d4 the option -Tlmit GE,LE will provide (Tgas>=10. AND Tgas<=1d4) as\
			the reaction range of validity. Operators opLow and opHigh must be one of the following: LE, GE, LT, GT.")
		self.parser.add_argument("-unsafe", action="store_true", help="skip to check if the build folder is empty or not")
		self.parser.add_argument("-useCustomCoe", help="use a user-defined custom function that returns a real*8 array of size\
			NREA = number of reactions, that replaces the standard rate coefficient calculation function. Note that FUNCTION\
			must be explicitly included in krome_user_commons module.", metavar="FUNCTION")
		self.parser.add_argument("-useDvodeF90", action="store_true", help="use Dvode implementation in F90 (slower)")
		self.parser.add_argument("-useEquilibrium", action="store_true", help="check if the solver has reached the equilbirum.\
			If so break the solver's loop and return the values found. It is useful when the system oscillates around\
			a solution (as in some photoheating cases). To be used with caution!")
		self.parser.add_argument("-useFileIdx", action="store_true", help="use the reaction index in the reaction file instead of\
			using the automatic progressive index starting from 1. Useful with rate coefficients that depends on other\
			coefficients, e.g. k(10) = 1d-2*k(3)")
		self.parser.add_argument("-useH2opacity", action="store_true",help="use H2 opacity for H2 cooling")
		self.parser.add_argument("-useIERR", action="store_true",help="use ierr in the interface with KROME to return errors instead\
			of stopping the exectution")
		self.parser.add_argument("-useN", action="store_true",help="use number densities (1/cm3) as input/ouput instead of fractions (#)")
		self.parser.add_argument("-useODEConstant", help="postpone an expression to each ODE. EXPRESSION must be a valid f90\
			expression (e.g. *3.d0 or +1.d-10)", metavar="EXPRESSION")
		self.parser.add_argument("-usePhIoniz", action="store_true", help="include photochemistry")
		self.parser.add_argument("-usePlainIsotopes", action="store_true", help="use kA format for isotopes instead of [k]A format,\
			where k is the isotopic number and A is the atom name, e.g. krome looks for 14C instead of [14]C in the reactions file.")
		self.parser.add_argument("-useTabs", action="store_true", help="use tabulated rate coefficients (free parameter: temperature)")
		self.parser.add_argument("-v", action="store_true", help="print the current version of KROME")
		self.parser.add_argument("-ver", action="store_true", help="same as -v")
		self.parser.add_argument("-version", action="store_true", help="same as -v")
	 
	
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
			[argv.append(x) for x in ["-useN","-iRHS","-skipJacobian"]]
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
			filename = "networks/react_primordial3"
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
		elif(args.test=="stars"):
			[argv.append(x) for x in ["-star","-usePlainIsotopes","-nomassCheck"]]
			filename = "networks/react_star"
		elif(args.test=="coolCO"):
			[argv.append(x) for x in [""]]
			filename = "networks/react_dummy"
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

		args = self.parser.parse_args() #return namespace from argv

		#use custom option file (load options from a file and append to argv)
		if(args.options):
			fopt = args.options.strip() #get filename
			print "Reading option -option="+fopt
			#check if option file exsists
			if(not(file_exists(fopt))):
				print "ERROR: custom option file \""+fopt+"\" does not exist!"
				sys.exit()
			#read from file
			fho = open(fopt,"rb")
			for row in fho:
				srow = row.strip()
				if(srow==""): continue #skip blank lines
				arow = srow.split()
				#append to argv
				for x in arow:
					sys.argv.append(x)

			args = self.parser.parse_args() #return updated namespace
		
		#get a citation and exit
		if(args.quote):
			get_quote()
			sys.exit()

		#save options into a file
		fopt = open("options.log","w")
		for k,v in vars(args).iteritems():
			#if option is set add to the namespace
			if(v):
				if(v is True): v="" #if is exactly True write key only
				fopt.write("-"+k+" "+v+"\n") #write to file

		#you can select only one -forceMF
		if((args.forceMF222) and (args.forceMF21)):
			die("ERROR: options -forceMF222 and -forceMF21 are mutually exclusive: choose one.")

		#get filename
		if(not(self.is_test) and args.n): self.filename = args.n
		#chech if reactions file exists
		if(args.n or self.is_test):
			if(not(os.path.isfile(self.filename))): die("ERROR: Reaction file \""+self.filename+"\" doesn't exist!")
		else:
			die("ERROR: you must define -n FILENAME, where FILENAME is the reaction file!")

		#read the coolFile
		if(args.coolFile):
			self.coolFile = args.coolFile
			print "Reading option -coolFile (filename="+str(self.coolFile)+")"


		#use f90 solver
		if(args.useDvodeF90):
			self.useDvodeF90 = True
			self.solver_MF = 227
			print "Reading option -useDvodeF90"

		#set implicit RHS
		if(args.iRHS):
			self.use_implicit_RHS = True
			self.solver_MF = 222
			if(self.useDvodeF90):
				self.solver_MF = 227
			print "Reading option -iRHS"
		#force MF=21
		if(args.forceMF21):
			self.solver_MF = 21
			if(self.useDvodeF90):
				self.solver_MF = 27
			print "Reading option -forceMF21"
		#force MF=222
		if(args.forceMF222):
			self.solver_MF = 222
			if(self.useDvodeF90):
				self.solver_MF = 227
			print "Reading option -forceMF222"
		#use numeric density instead of fractions as input
		if(args.useN):
			self.useX = False
			print "Reading option -useN"

		#use rate tables
		if(args.useTabs):
			self.useTabs = True
			print "Reading option -useTabs"

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

                #determine H2shielding types 
		if(args.shielding):
			myShielding = [x.strip() for x in args.shielding.split(",")]
			#list of the shielding approximations 
			allShielding = ["DB96","WG11"]
			for shi in myShielding:
				if(not(shi in allShielding)):
					die("ERROR: Shielding \""+shi+"\" is unknown!\nAvailable shielding are: "+(", ".join(allShielding)))
			if(len(myShielding)>1): die("ERROR: "+(", ".join(allShielding))+" are mutually exclusive!")
			self.useShieldingDB96 = ("DB96" in myShielding)
			self.useShieldingWG11 = ("WG11" in myShielding)
                        self.useShielding = True
			print "Reading option -shielding (TYPE="+(",".join(myShielding))+")"


		#use human Fluxes
		if(args.compressFluxes):
			self.humanFlux = False
			print "Reading option -compressFluxes"

		#use cooling dT/dt in the ODE fex
		if(args.skipODEthermo):
			self.useODEthermo = False
			print "Reading option -skipODEthermo"

		#use species mass conservation (and charge)
		if(args.conserve):
			self.useConserve = True
			self.useConserveE = True
			print "Reading option -conserve"

		#use species charge conservation only
		if(args.conserveE):
			self.useConserveE = True
			print "Reading option -conserveE"

		#same index for equivalent reactions with different Tlimits
		if(args.mergeTlimits):
			self.mergeTlimits = True
			print "Reading option -mergeTlimits"

		#use short header for f90 files
		if(args.sh):
			self.shortHead = True
			print "Reading option -sh"

		#enable thermochemistry checking
		if(args.checkThermochem):
			self.checkThermochem = True
			print "Reading option -checkThermochem"

		#use IERR interface for krome
		if(args.useIERR or args.ierr):
			self.useIERR = True
			print "Reading option -useIERR"

		#check if reverse reactions are present in the network
		if(args.checkReverse):
			self.checkReverse = True
			print "Reading option -checkReverse"

		#do not write anything to the build directory
		if(args.dry):
			self.isdry = True
			print "Reading option -dry"

		#use short header for f90 files
		if(args.v or args.ver or args.version):
			print "You are using KROME "+self.version+" \""+self.codename+"\""
			#print "MD5: "+GetHashofDirs()
			sys.exit()

		#skip reaction mass / charge check
		if((args.nomassCheck and args.nochargeCheck) or args.noCheck):
			print "Reading option -nochargeCheck"
			print "Reading option -nomassCheck"
			self.checkMode = "NONE"
		elif(args.nomassCheck and not(args.nochargeCheck)):
			print "Reading option -nomassCheck"
			self.checkMode = "CHARGE"
		elif(not(args.nomassCheck) and (args.nochargeCheck)):
			print "Reading option -nochargeCheck"
			self.checkMode = "MASS"
		elif(not(args.nomassCheck) and not(args.nochargeCheck)):
			self.checkMode = "ALL"
		else:
			print "ERROR: problem with -nomassCheck and/or -nochargeCheck and/or -noCheck"
			sys.exit()

		#use nuclear multeplicity flux/(1.+delta_ij)
		if(args.nuclearMult):
			self.useNuclearMult = True
			print "Reading option -useNuclearMult"

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
			if(args.heating):
				if("COMPR" in args.heating):
					print "ERROR: -heating=COMPRESS is intended only for one-zone gravitational collapse!"
					sys.exit()

			if(not(args.customATOL) and not(args.ATOL)):
				print "WARNING: default ATOL set to 1e-10 due to -ramses flag."
				self.ATOL = 1e-10

		#creates ramsesTH version patches
		if(args.ramsesTH):
			self.doRamsesTH = True
			print "Reading option -ramsesTH"
			if(not(args.compact)):
				die("ERROR: the patch for RAMSES TH requires the -compact option!")
			if(self.is_test):
				die("ERROR: -test option and -ramsesTH are incompatible!")
			if(self.useX):
				die("ERROR: the patch for RAMSES TH requires the -useN option!")
			if(args.heating):
				if("COMPR" in args.heating):
					die("ERROR: -heating=COMPRESS is intended only for one-zone gravitational collapse! Remove it")

			if(not(args.customATOL) and not(args.ATOL)):
				print "WARNING: default ATOL set to 1e-10 due to -ramsesTH flag."
				self.ATOL = 1e-10

		#creates flash patches
		if(args.flash):
			self.doFlash = True
			print "Reading option -flash"
			if(not(args.compact)):
				print "ERROR: the patch for FLASH requires the -compact option!"
				sys.exit()
			if(self.is_test):
				print "ERROR: -test option and -flash are incompatible!"
				sys.exit()
			if(self.useX):
				print "ERROR: the patch for FLASH requires the -useN option!"
				sys.exit()
			if(args.heating):
				if("COMPR" in args.heating):
					print "ERROR: -heating=COMPRESS is intended only for one-zone gravitational collapse!"
					sys.exit()

			if(not(args.customATOL) and not(args.ATOL)):
				print "WARNING: default ATOL set to 1e-10 due to -flash flag."
				self.ATOL = 1e-10


		#creates enzo patches
		if(args.enzo):
			self.doEnzo = True
			print "Reading option -enzo"
			if(not(args.compact)):
				print "ERROR: the patch for ENZO requires the -compact option!"
				sys.exit()
			if(self.is_test):
				print "ERROR: -test option and -enzo are incompatible!"
				sys.exit()
			if(self.useX):
				print "ERROR: the patch for ENZO -useN option!"
				sys.exit()
			if(args.heating):
				if("COMPR" in args.heating):
					print "ERROR: -heating=COMPRESS is intended only for one-zone gravitational collapse!"
					sys.exit()
			if(not(args.customATOL) and not(args.ATOL)):
				print "WARNING: default ATOL set to 1e-10 due to -enzo flag."
				self.ATOL = 1e-10

		#creates simple C wrapper
		if(args.C):
			self.wrapC = True
			print "Reading option -C"


		#skip writing Jacobian in krome_ode.f90, allows faster compilation
		if(args.skipJacobian):
			self.doJacobian = False
			print "Reading option -skipJacobian"


		#skip checking for objects in build/
		if(args.unsafe):
			self.safe = False
			print "Reading option -unsafe"

		#enable stellar physics
		if(args.stars):
			self.useStars = True
			print "Reading option -stars"

		#do not write test.f90 or Makefile
		if(args.noExample):
			self.noExample = True
			print "Reading option -noExample"


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
			myCools = [x.strip() for x in myCools]
			#list of all cooling (excluded from file)
			allCools = ["ATOMIC","H2","HD","DH","DUST","H2GP98","COMPTON","EXPANSION","CIE","CONT","CHEM","DISS","Z"]
			fileCools = [] #list of the cooling read from file
			#load additional coolings from file
			fname = self.coolFile
			if(not(file_exists(fname))):
				print "ERROR: file "+fname+" not found!"
				sys.exit()
			fh = open(fname,"rb")
			inComment = False
			for row in fh:
				srow = row.strip()
				#skip cooments
				if(srow==""): continue
				if(srow[0]=="#"): continue
				if(srow[:1]=="//"): continue
				if(srow[:1]=="/*"): inComment = True
				if("*/" in srow): 
					inComment = False
					continue
				srow = srow.split("#")[0]
				if(inComment): continue
				#look for the metal name
				if("metal:" in srow):
					metal_name = srow.split(":")[1].strip().capitalize() #metal name starts with capital letter
					#metal name in roman style (e.j. C+ -> CII)
					if("+" in metal_name):
						mname = metal_name.replace("+","") + int_to_roman(metal_name.count("+")+1)
					elif("-" in metal_name):
						mname = metal_name.replace("-","") + "m"+int_to_roman(metal_name.count("-")+1)
					else:
						mname = metal_name+"I" #is not an ion

					if(mname in allCools):
						print "ERROR: conflict name for "+mname+", which is already present!"
						sys.exit()
					fileCools.append({"flag":mname,"name":metal_name}) #append to the list fo the available coolings
					allCools.append(mname) #append flag to the list of the coolants
			#write found coolants
			if(len(fileCools)>0): print "Cooling "+(", ".join([x["flag"] for x in fileCools]))+" available from "+fname
			
			#check coolant names
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
			if("EXPANSION" in myCools): self.useCoolingExpansion = True
			if("CHEM" in myCools): self.useCoolingChem = True
			if("CIE" in myCools): self.useCoolingCIE = True
			if("DISS" in myCools): self.useCoolingDISS = True
			if("CONT" in myCools): self.useCoolingCont = True
			if("Z" in myCools): self.useCoolingZ = True

			#loop over metals loaded from file and search for themin the cooling flags provided by the user
			for met in fileCools:
				if(met["flag"] in myCools):
					#print "Option "+met["flag"]+" will load data from "+fname
					self.useCoolingZ = True
					self.Zcools.append(met["name"]) #append metal name to the list of the requested species
					self.zcoolants.append(met["flag"]) #append roman metal name to the list of the coolants
					#add also neutral in case of ions
					if("+" in met["name"]):
						neutral_name = met["name"].replace("+","")
						if(not(neutral_name in self.Zcools)): self.Zcools.append(neutral_name)
					#add also neutral in case of anions
					if("-" in met["name"]):
						neutral_name = met["name"].replace("-","")
						if(not(neutral_name in self.Zcools)): self.Zcools.append(netural_name)

			self.use_cooling = True
			self.hasDust = False
			for aa in argv:
				if("dust=" in aa): self.hasDust = True
			if(self.useCoolingDust and not(self.hasDust)):
				die("ERROR: to include dust cooling you need dust (use -dust=[see help]).")
			if(("CHEM" in myCools) and ("ATOMIC" in myCools)):
				die("ERROR: CHEM and ATOMIC cooling are mutually exclusive!")
			if(("CIE" in myCools) and ("CONT" in myCools)):
				die("ERROR: CIE and CONT cooling are mutually exclusive!")

			self.use_thermo = True

			print "Reading option -cooling ("+(",".join(myCools))+")"
		
		#determine heating types
		if(args.heating):
			myHeat = args.heating.upper().split(",")
			myHeat = [x.strip() for x in myHeat]
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
		if(args.forceRWORK):
			myrwork = args.forceRWORK
			self.force_rwork = True
			print "Reading option -forceRWORK (RWORK="+str(myrwork)+")"

		#use custom function for coefficient instead of coe_tab()
		if(args.useCustomCoe):
			self.customCoeFunction = args.useCustomCoe
			self.useCustomCoe = True
			print "Reading option -useCustomCoe (Expression="+str(self.customCoeFunction)+")"

		#use function to append after each ODE
		if(args.useODEConstant):
			self.ODEConstant = args.useODEConstant
			self.useODEConstant = True
			print "Reading option -useODEConstant (Constant="+str(self.ODEConstant)+")"

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

		#project name folder
		if(args.source):
			flist = ["krome_commons.f90", "krome_cooling.f90", "krome.f90", "krome_heating.f90"]
			flist += ["krome_photo.f90","krome_subs.f90", "krome_user_commons.f90", "krome_constants.f90"]
			flist += ["krome_dust.f90", "kromeF90.f90", "krome_ode.f90", "krome_reduction.f90"]
			flist += ["krome_tabs.f90", "krome_user.f90"]

			src = str(args.source)
			print "Reading option -source (name="+src+")"
	
			#check if folder exists
			if(not(os.path.exists(src))):
				print "ERROR: the folder "+src+"/ doesn't exist!"
				sys.exit()
			#check if the file in flist are present in the folder
			notfound = []
			for fle in flist:
				if(not(file_exists(src+"/"+fle))):
					notfound.append(fle)

			#if file missing write the error message
			if(len(notfound)>0):
				print "ERROR: you suggested to use the folder "+src+"/ as source"
				print " but the following file(s) missing:"
				print " " + (", ".join(notfound))
				sys.exit()

			self.srcFolder = src+"/"


		#typeGamma
		if(args.gamma):
			typeGamma = args.gamma
			self.typeGamma = typeGamma.replace("\"","")
			if(not(args.heating) and not(args.cooling)):
				print "ERROR: you are trying to use -gamma without -cooling or -heating"
				sys.exit()
			print "Reading option -gamma (gamma="+str(self.typeGamma)+")"



		#offset for ramses
		if(args.ramsesOffset):
			if(not(args.ramses)): die("ERROR: if you use -ramsesOffset you must also add -ramses option!")
			self.ramses_offset = args.ramsesOffset
			print "Reading option -ramsesOffset (offset="+str(args.ramsesOffset)+")"

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

	###################################################
	def safe_check(self):
		if(not(self.safe) or not(os.path.exists(self.buildFolder))): return
		if(self.isdry): return
		wlk = os.walk(self.buildFolder).next()
		wlk = wlk[1]+wlk[2] #folders+files
		if(len(wlk)<1): return
		print "************************************************"
		print "WARNING: the folder "+self.buildFolder+" is not empty"
		print " some items may be replaced. Do you want to proceed?"
		print "To avoid this message use -unsafe option."
		print "************************************************"
		a = raw_input("Any key to ignore q to quit... ")
		if(a=="q"): print sys.exit()

	####################################################
	#load thermochemistry data from chemkin-formatted file
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

		#exited levels of some molecules add here if needed
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
		TminAuto = self.TminAuto
		TmaxAuto = self.TmaxAuto
		hasFormat = False
		format_items = 4+len(ireact)+len(iprod)
		if(skipDup): fdup = open("duplicates.log","w")
		idxFound = tminFound = tmaxFound = rateFound = True
		qeffFound = False
		group = "__DEFAULT__" #default group for reactions
		specs = []
		reacts = []
		reags = [] #list of reagents for already found
		prods = [] #list of prods for already found
		idxs = [] #list of index for already found

		#read the size of the file in lines (skip blank and comments)
		# to have a rough idea of the size
		fh = open(filename,"rb")
		line_count = 0
		for row in fh:
			if(row.strip()==""): continue
			if(row.strip()[0]=="#"): continue
			line_count += 1
		fh.close()
	
		#warning if the number of lines exceed a certain limit
		if(line_count>1000): print "Found "+str(line_count)+" lines! It takes a while..."

		#star proper reading
		fh = open(filename,"rb") #OPEN FILE
		isComment = False #flag for comment block
		noTabNext = False #flag for use tabs for the next reaction 
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

			#search for group indication
			if("@group:" in srow):
				group = srow.replace("group:","").strip().replace(" ","_")
				if(not(group.isalnum())): die("ERROR: group must be alphanumeric. Found "+group)
				print "Found reactions group "+group
				if(not(group in self.groups)): groups.append(group)
				continue

			#search for variables
			if("@var:" in srow):
				arow = srow.replace("@var:","").split("=")
				if(len(arow)!=2):
					print "ERROR: variable line must be @var:variable=F90_expression"
					print "found: "+srow
					sys.exit()
				if(arow[0] in self.coevars): continue #skip already found variables
				self.coevars[arow[0]] = [ivarcoe,arow[1]]
				ivarcoe += 1 #count variables to sort
				continue #SKIP: a variable line is not a reaction line

			#search for ghost species
			if("@ghost:" in srow):
				ghost = srow.replace("@ghost:","").strip()
				print "Found ghost species: "+ghost
				mol = parser(ghost,mass_dic,atoms,self.thermodata)
				if(not(mol.name in spec_names)):
					spec_names.append(mol.name)
					specs.append(mol)
				mol.idx = spec_names.index(mol.name) + 1
				continue #SKIP: a ghost line is not a reaction line
				

			#search for format string
			if("@format:" in srow):
				idxFound = tminFound = tmaxFound = rateFound = qeffFound = False
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
					if(x=="qeff" or x=="qpp"):
						iqeff = ipos
						qeffFound = True
					ipos += 1 #increase position
				#check for rate
				if(not(rateFound)):
					print "ERROR: format must contain rate token"
					sys.exit()
			
				continue #SKIP format line (it is not a reaction line)

			#if requested the next reaction will not uses tabs
			if(srow.lower()=="@notabnext"):
				noTabNext = True
				continue #SKIP (not a reaction)
				

			arow = srow.split(self.separator,format_items-1) #split only N+1 elements with N seprations
			arow = [x.strip() for x in arow] #strip single elements
			if(len(arow)!=format_items): 
				print "WARNING: wrong format for reaction "+str(rcount+1)
				print srow
				a = raw_input("Any key to continue q to quit... ")
				if(a=="q"): print sys.exit()
				continue #check line format (N elements, 4=idx+Tmin+Tmax+rate)
			found_one = True #flag to determine at least one reaction found
			rcount += 1 #count the totoal number of reaction found
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

			reactants = [arow[x].strip().upper() for x in ireact]
			products = [arow[x].strip().upper() for x in iprod]

			#store reactants "curlyness" before purge
			myrea.curlyR = [("{" in x) for x in reactants]

			#purge reactants name from curly brackets
			reactants = [x.replace("}","").replace("{","") for x in reactants]

			#remove empty reactants/products (i.e. dummy) and sort
			reags_clean = sorted([x for x in reactants if x.strip()!=""])
			prods_clean = sorted([x for x in products if x.strip()!=""])

			#check for identical ractions
			foundAlready = False
			for i in range(len(reags)):
				if(reags_clean==reags[i] and prods_clean==prods[i] and self.mergeTlimits):
					foundAlready = True
					rcount -= 1 #decrease reaction index (since already increased few lines above)
					myrea.idx = idxs[i]
					print "already found: ("+str(idxs[i])+") "+(" + ".join(reags_clean))+" -> "+(" + ".join(prods_clean))

			#store reactants and products to find identical reactions (and different Tlimits)
			if(not(foundAlready)):
				reags.append(reags_clean)
				prods.append(prods_clean)
				idxs.append(myrea.idx)

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

			myrea.Tmin = "2.73d0" #default min temperature
			myrea.Tmax = "1.d8" #default max temperature
			#search for reactions without Tlims
			if(tminFound):
				if(arow[iTmin].strip().upper() in ["N","NONE","N/A","NO",""]): myrea.hasTlimitMin = False
			if(tmaxFound):
				if(arow[iTmax].strip().upper() in ["N","NONE","N/A","NO",""]): myrea.hasTlimitMax = False
			#store Tlimits if any
			if(myrea.hasTlimitMin):
				if(tminFound): myrea.Tmin = format_double(arow[iTmin]) #get Tmin
				if(tminFound): TminAuto = min(float(arow[iTmin].lower().replace("d","e")), TminAuto)
			if(myrea.hasTlimitMax):
				if(tmaxFound): myrea.Tmax = format_double(arow[iTmax]) #get Tmax
				if(tmaxFound): TmaxAuto = max(float(arow[iTmax].lower().replace("d","e")), TmaxAuto)
			#store other data
			area = arow[irate].split(":",2) #raction can be if_condition:reaction_rate
			if(len(area)==1):
				myrea.ifrate = "" #store empty prepending if condition
				myrea.krate = arow[irate] #get reaction rate written in F90 style
			else:
				myrea.ifrate = area[0] #store prepending if condition
				myrea.krate = area[1] #get reaction rate written in F90 style
				
			if(qeffFound): myrea.qeff = arow[iqeff]

			#if(self.useCustomCoe): myrea.krate = "0.d0" #when custom function is used standard coefficient are set to zero
			#loop over reactants to grep molecules
			for r in reactants:
				if(r.strip()=="G" and not(self.use_photons)): continue
				if(r.strip()=="E-"): r = "E"
				if(r.strip()!=""):
					mol = parser(r,mass_dic,atoms,self.thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
					myrea.reactants.append(mol) #add molecule object to reactants
			#loop over prodcuts to grep molecules
			for p in products:
				if(p.strip()=="G" and not(self.use_photons)): continue
				if(p.strip()=="E-"): p = "E"
				if(p.strip()!=""):
					mol = parser(p,mass_dic,atoms,self.thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
					myrea.products.append(mol) #add molecule object to products

			myrea.build_verbatim() #build reaction as string (e.g. A+B->C)
			#myrea.reactants = sorted(myrea.reactants, key=lambda r:r.idx) #sort reactants
			#myrea.products = sorted(myrea.products, key=lambda p:p.idx) #sort products
			myrea.build_RHS(self.useNuclearMult) #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
			myrea.build_phrate() #build photoionization rate
			myrea.check(self.checkMode) #check mass and charge conservation
			myrea.group = group #add the group to the reaction
			myrea.canUseTabs = not(noTabNext) #check if this reaction can use tabs or not
			if(myrea.krate.count("(")!=myrea.krate.count(")")):
				print "ERROR: unbalanced brakets in reaction "+str(myrea.idx)
				print " "+myrea.verbatim
				print " rate = "+myrea.krate
				sys.exit()

			#skip duplicated reactions if requested
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

			#append reactions if not skipped 
			if(not(skip_append)): reacts.append(myrea)
			noTabNext = False #return to default value

	

		if(skipDup): 
			fdup.close()
			print "Skipped duplicated reactions:",skipped_dupl

		#check file format
		if(not(found_one)):
			die("ERROR: no valid reactions found in file \""+filename+"\"")
		if(unmatch_idx):
			print "WARNING: index in \""+filename+"\" are not sequential!"

		#count reactions with unique index
		idxs = []
		nrea = 0
		for rea in reacts:
			if(rea.idx in idxs): continue #skip reactions same index
			idxs.append(rea.idx)
			nrea += 1

		#copy local to global vars
		self.nrea = nrea
		self.specs = specs
		self.reacts = reacts
		self.TminAuto = TminAuto
		self.TmaxAuto = TmaxAuto

		print "done!"
	
	#####################################################
	def photo_warnings(self):
		if(self.is_test): return #skip warning if test mode
		if(self.usePhIoniz or self.useHeatingPhoto):
			print "************************************************"
			print "REMINDER: note that, since you are using photon-based"
			print " options, you need to include \"call krome_init_photo()\""
			print " to your main file in order to initialize the rate" 
			print " coefficients!"
			print "************************************************"
			a = raw_input("Any key to continue...")
		
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
					myrev.hasTlimtMin = myrea.hasTlimitMin #get info on TlimitMin
					myrev.hasTlimtMax = myrea.hasTlimitMax #get info on TlimitMax
					myrev.reactants = myrea.products
					myrev.products = myrea.reactants
					myrev.TminOp = myrea.TminOp #get Tmin operator
					myrev.TmaxOp = myrea.TmaxOp #get Tmax operator
					myrev.build_verbatim() #build reaction as string (e.g. A+B->C)
					myrev.verbatim += " (REV)" #append REV label to reaction string
					myreactants = myrev.reactants #get list of reactants
					myproducts = myrev.products #get list of products
					myrev.curlyR = [False]*len(myrev.reactants)
					myrev.krate = myrea.krate #set the forward reaction
					myrev.krate = myrev.doReverse() #compute reverse reaction using Simoncini2013PRL
					myrev.build_RHS(self.useNuclearMult) #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
					myrev.check(self.checkMode) #check mass and charge conservation
					reacts.append(myrev)
			print "Inverse reaction added: "+str(count_reverse)
			
			self.reacts = reacts
			self.nrea = len(reacts)

	#########################
	#break degenerancy for reactions that contains catalyser
	# e.g. H + e -> H+ + e + e, which is equal to H -> H+ + e
	# arg1 and arg2 are the arrays of reactants and products
	# to return pruned reactant list
	def PRuniq(self,arg1, arg2):
		
		arg1n = []
		arg1c = [x.name for x in arg1]
		arg2c = [x.name for x in arg2]
		for x in arg1c:
			if(x in arg2c):
				ii = arg2c.index(x)
				arg2c[ii] = "@@@@"
				continue
			arg1n.append(x)
		return arg1n
			
	###########################################
	#check if reactions have their reverse in the chemical network
	def check_reverse(self):
		if(not(self.checkReverse)): return
		idxRev = []
		reacts = self.reacts
		#loop on reacts
		for i in range(len(reacts)):
			if(i in idxRev): continue #if already found reverse skip
			rea1 = reacts[i]
			#prune catalysers from reactants and products
			R1 = self.PRuniq(rea1.reactants, rea1.products)
			P1 = self.PRuniq(rea1.products, rea1.reactants)
			revFound = False
			#loop on reacts
			for j in range(i,len(reacts)):
				rea2 = reacts[j]
				#prune catalysers from reactants and products
				R2 = self.PRuniq(rea2.reactants, rea2.products)
				P2 = self.PRuniq(rea2.products, rea2.reactants)
				#verify conditions
				if((R1==P2) and (R2==P1)):
					revFound = True
					#store reverse index
					idxRev.append(i)
					idxRev.append(j)
					break #break loop when reverse found
			if(not(revFound)): print "WARNING: no reverse reaction found for "+rea1.verbatim
					
	###################################################
	def verifyThermochem(self):
		if(not(self.checkThermochem)): return
		for x in self.specs:
			if(not(x.name in self.thermodata)):
				print "WARNING: no thermochemical data for "+x.name+"!"



	###########################################add all the metals to cooling
	def addMetals(self):
		Zcools = self.Zcools
		specs = self.specs
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
				print "Adding species \""+zcool+"\" (requested by metal cooling)"
				mymol = parser(zcool,self.mass_dic,self.atoms,self.thermodata)
				mymol.idx = len(specs)+1
				self.specs.append(mymol)
		self.totMetals = "tot_metals = " + (" + ".join(["n(idx_"+x.replace("+","j")+")" for x in Zcools]))
	
	#######################################
	def addReaMin(self):
		for rea in self.reacts:
			if(not("krome_" in rea.krate)): rea.krate = "small + ("+rea.krate+")"


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
		#	if(mol.name=="G"):
		#		has_g = True
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
		fout.write("extern int krome_nrea = "+str(self.nrea)+"; //number of reactions\n")
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
		print "Species list saved in species.log"
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
		idxs = [] #already employed indexes
		for rea in reacts:
			if(rea.idx in idxs): continue #skip if already employed index
			idxs.append(rea.idx)
			rhs = "kflux("+str(rea.idx)+")"
			if(self.humanFlux): rhs = rea.RHS
			for r in rea.reactants:
				dns[r.idx-1] = dns[r.idx-1].replace(" = 0.d0"," =")
				dns[r.idx-1] += " -"+rhs
			for p in rea.products:
				dns[p.idx-1] = dns[p.idx-1].replace(" = 0.d0"," =")
				dns[p.idx-1] += " +"+rhs
				


		#equilibrium matrix (NOTE: NOT SUPPORTED!)
		#TODO: fix it
		if(False):
			idxs = [] #already employed indexes
			xvar = [] #names of the composite variables
			xvarM = [] #equilibrium matrix
			xvarS = [] #symbolic matrix (for signless comparison)
			for rea in reacts:
				if(rea.idx in idxs): continue #skip if already employed index
				ridx = rea.idx-1
				#build matrix variable reactants
				pre_xvarR = []
				for r in rea.reactants:
					pre_xvarR.append("n("+r.fidx+")")
				#build matrix variable reactants
				pre_xvarP = []
				for p in rea.products:
					pre_xvarP.append("n("+p.fidx+")")
				#look for powers (e.g. H*H*H, or H2*H2)
				#equals = 0
				#if(pre_xvar.count(pre_xvar[0])==len(pre_xvar)): equals = pre_xvar.count(pre_xvar[0])
				pre_xvarR = "*".join(pre_xvarR)
				pre_xvarP = "*".join(pre_xvarP)
				#in case of new variable add a row to the matrix
				if(not(pre_xvarR in xvar)):
					xvar.append(pre_xvarR)
					xvarM.append(["0d0" for i in reacts])
					xvarS.append(["" for i in reacts])
				if(not(pre_xvarP in xvar)):
					xvar.append(pre_xvarP)
					xvarM.append(["0d0" for i in reacts])
					xvarS.append(["" for i in reacts])
				#add rate to the matrix element
				for r in rea.reactants:
					xvarM[xvar.index(pre_xvarR)][ridx] = xvarM[xvar.index(pre_xvarR)][ridx].replace("0d0", "")
					xvarM[xvar.index(pre_xvarR)][ridx] += " -k("+str(ridx+1)+")"
					xvarS[xvar.index(pre_xvarR)][ridx] += "_"+str(ridx)
				for p in rea.products:
					xvarM[xvar.index(pre_xvarP)][ridx] = xvarM[xvar.index(pre_xvarP)][ridx].replace("0d0", "")
					xvarM[xvar.index(pre_xvarP)][ridx] += " +k("+str(ridx+1)+")"
					xvarS[xvar.index(pre_xvarP)][ridx] += "_"+str(ridx)

			#add conseravtion to matrix
			xvarM.append(["1d0" for i in reacts])
			xvarS.append(["eq1" for i in reacts]) #add dummy for symbolic matrix
			xvar.append("+".join(xvar))
		
			#remove linear dependent terms
			xvarM_u = []
			for i in range(len(xvarS)):
				row1 = xvarS[i]
				isequal = False
				for j in range(i+1,len(xvarS)):
					row2 = xvarS[j]
					if(row1==row2):
						isequal = True
						break
				if(not(isequal)): 
					xvarM_u.append(xvarM[i])
			if(len(xvar)==len(xvarM_u)): print "NOTE: this system can be solved algebrically to the equilibrium"
			#print str(len(xvar))+" variables and "+str(len(xvarM_u))+ " equations"

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
		idxs = []
		for rea in reacts:
			if(rea.idx in idxs): continue #avoid reactions with same index
			idxs.append(rea.idx)
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
		idxs = [] #store reaction indexes
		for rea in reacts:
			if(rea.idx in idxs): continue #skip reactions with same index
			idxs.append(rea.idx) #store reaction index
			#loop over reactants
			for ri in range(len(rea.reactants)):
				if(rea.curlyR[ri] and self.useNuclearMult): continue #skip curly reactants
				sjac = rea.nuclearMult+"k("+str(rea.idx)+")" #init and include nuclearMulteplicity if any
				r1 = rea.reactants[ri] #get ri-th reactant
				#loop over reactants again
				for rj in range(len(rea.reactants)):
					if(rea.curlyR[rj] and self.useNuclearMult): continue #skip curly reactants
					r2 = rea.reactants[rj] #get rj-th reactant
					#if reactants are different add rj-th to jacobian
					if(ri!=rj):
						sjac += "*n("+str(r2.fidx)+")"
				#update built jacobian and sparsity for reactants
				for rr in rea.reactants:
					jac[rr.idx-1][r1.idx-1] = jac[rr.idx-1][r1.idx-1].replace(" = 0.d0"," =")
					jac[rr.idx-1][r1.idx-1] += " &\n-"+sjac
					jsparse[rr.idx-1][r1.idx-1] = 1
				#update built jacobian and sparsity for products
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
	#this function computes the sparsity structure in the yale format
	# namely IA e JA in the DLSODES documentation
	#Do the same for dvodeF90
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

		#sparsity for dvodeF90
		jaca = [] #unrolled sparse jacobian
		if(self.useDvodeF90):
			ia = [1] #sparsity structure IA (see dvode_f90 manual)
			ja = [] #sparsity structure JA (see dvode_f90 manual)
			isum = 1
			for i in range(neq):
				isum += sum(jsparse[i])
				ia.append(isum)
				for j in range(neq):
					if(jsparse[i][j]==1):
						ja.append(j+1)
						jaca.append(self.jac[i][j].replace("pdj("+str(i+1)+") =  &\n","").strip())
			iaf = "iauser(:) = (/" + (",".join([str(x) for x in ia])) + "/)"
			jaf = "jauser(:) = (/"+(",".join([str(x) for x in ja]))+"/)"


		#if MF=222 no need for sparsity structure arrays
		if(self.solver_MF == 222):
			iaf = jaf = ""
		self.ja = ja
		self.ia = ia
		self.iaf = iaf
		self.jaf = jaf
		self.jaca = jaca

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
		nrea = self.nrea #number of reactions
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
		elif(solver_miter==7):
			lwm = 0
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
		elif(solver_MF==227):
			lrw = 0
		elif(solver_MF==27):
			lrw = 0
		else:
			die("ERROR: solver_MF value "+str(solver_MF)+" unknown in LRW calculation!")

		if(self.force_rwork):
			lrw = myrwork
		#lrw = int(20+(2+0.5)*nnz + (11+4.5)*neq) #RWORK size

		print " LWM:",lwm,"LRW:",lrw
		self.lwm = lwm
		self.lrw = lrw

	#######################################
	#this function loads the cooling functions from a file
	# and prepares all the necessary stuff 
	def createZcooling(self):
		fname = self.coolFile #file containing cooling information
		if(not(file_exists(fname))):
			print "ERROR: file "+fname+" not found!"
			sys.exit()
	
		zcoolants = self.zcoolants
		skip = False
		inmetal = False
		skip_metal = False
		vars_cool = [] #variables
		coolZ_functions = [] #cooling functions
		krates = [] #rates
		nkrates = 0
		needOrthoPara = False
		#read file
		print "******************"
		print "Reading coolants from "+fname+"..."
		fh = open(fname,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue #skip blank lines
			if(srow[0]=="#"): continue #skip lines with coments
			if(srow[:1]=="//"): continue #comments
			if(srow[:1]=="/*"): skip = True #skip large comments
			if("*/" in srow): 
				skip = False
				continue
			srow = srow.split("#")[0] #skip comments
			if(skip): continue
			if((srow=="endmetal" or srow=="end metal") and skip_metal):
				inmetal = False
				skip_metal = False
				continue
			if(skip_metal): continue #skip metals that are not present in the cooling options
			#read variables and store to vars_cool
			if("@var:" in srow):
				srow = srow.replace("@var:","").strip()
				vars_cool.append([x.strip() for x in srow.split("=")]) #read extra variables
				continue
			#read metal name
			if("metal:" in srow):
				inmetal = True #flag reading metal
				real_variables = []
				kkrates = []
				#modify the name of the metal (e.g. C+ becomes CII)
				metal_name = srow.replace("metal:","").strip()
				if("+" in metal_name):
					mname = metal_name.replace("+","") + int_to_roman(metal_name.count("+")+1)
				elif("-" in metal_name):
					mname = metal_name.replace("-","") + "m"+int_to_roman(metal_name.count("-")+1)
				else:
					mname = metal_name+"I" #is not an ion
				if(not(mname.upper() in [x.upper() for x in zcoolants])): skip_metal = True

				function_name = "cooling"+mname #name of the cooling function
				cur_metal = metal_name.replace("+","j").replace("-","k") #current metal name
				excitation_rates = [] #list of excitation rates written in F90 for cur_metal 
				transitions = [] #list of transitions
				colliders = [] #list of colliders names
				Aijs = dict() #init dictionary for Einsten's
				levels = dict() #init dictionary for levels
				#prepares the header of the function
				full_cool = "\n"
				full_cool += "!########## " + cur_metal + " #########\n"
				full_cool += "function "+function_name+"(n,inTgas,k)\n"
				full_cool += "use krome_commons\n"
				full_cool += "implicit none\n"
				full_cool += "real*8::"+function_name+",n(:),inTgas,Tgas,k(:),invTgas\n"
				full_cool += "#KROME_replace_with_declarations\n\n"
				full_cool += "Tgas = inTgas\n"
				full_cool += "invTgas = 1d0/Tgas\n"

			#read the levels
			if("level" in srow):
				nlev = int(srow.split(":")[0].replace("level","").strip()) #level number
				energy, gmult = srow.split(":")[1].split(",") #get energy and degenerancy factor
				levels[nlev] = {"energy":float(energy), "gmult":float(gmult)} #store energy and degenerancy
				continue
			
			#end of metals block
			if(srow=="endmetal" or srow=="end metal"):
				#check end metal
				if(not(inmetal)):
					print "ERROR: end metal statement without metal: statement!"
					sys.exit()
				inmetal = False #in metal block flag

				#prepare ortho/para variables if needed
				if(needOrthoPara):
					full_cool += "H2or = "+str(ortho2para/(1e0+ortho2para)) + " * n(idx_H2)\n"
					full_cool += "H2pa = "+str(1e0/(1e0+ortho2para)) + " * n(idx_H2)\n"
					full_cool += "!explicit if needed\n"
					full_cool += "!H2or = o2p / (1d0 + o2p) * n(idx_H2)\n"
					full_cool += "!H2pa = 1d0 / (1d0 + o2p) * n(idx_H2)\n"
			
				#write de-excitation rates
				full_cool += "\n!de-excitation rates\n"
				for x in kkrates:
					full_cool += x + "\n"
				full_cool += "\n"

				#prepares conversion factor to ji from ij
				full_cool += "!excitation rates\n"
				trfound = []
				for tr in transitions:
					ij2jivar = "g"+cur_metal+"_g"+str(tr["up"])+str(tr["down"])+"to"+str(tr["down"])+str(tr["up"]) 
					ij2ji = ij2jivar + " = " + str(float(levels[tr["up"]]["gmult"]) / levels[tr["down"]]["gmult"]) + "d0"
					ij2ji += " * exp(-" +str(float(levels[tr["up"]]["energy"]) - levels[tr["down"]]["energy"]) + "*invTgas)"
					#skip already found transitions					
					if(not([tr["up"], tr["down"]] in trfound)):
						full_cool += ij2ji + "\n"
						real_variables.append(ij2jivar)
						trfound.append([tr["up"], tr["down"]])
				
				#prepares rate conversion using conversion factor computed above
				full_cool += "\n"
				for tr in transitions:
					coll = tr["coll"]
					ij2jivar = "g"+cur_metal+"_g"+str(tr["up"])+str(tr["down"])+"to"+str(tr["down"])+str(tr["up"]) 
					varup = "g"+str(tr["up"])+str(tr["down"])+cur_metal+"_"+coll
					vardown = "g"+str(tr["down"])+str(tr["up"])+cur_metal+"_"+coll
					real_variables.append(varup)
					real_variables.append(vardown)
					full_cool += vardown + " = " + varup +" * " + ij2jivar + "\n"

				#prepares the Mij and Mji transition elements
				full_cool += "\n"
				full_cool += "!transitions\n"
				varMexist = []
				MM = dict()
				MMij = dict()
				#build Mij
				for tr in transitions:
					varM = "M" + str(tr["down"]) + str(tr["up"]) + cur_metal #init variable name
					real_variables.append(varM) #add to double variable list
					#check if variable is unique
					if(not(varM in varMexist)):
						MM[varM] = []
						MMij[varM] = [tr["down"],tr["up"]]
					coll = tr["coll"] #get colliders for the given transition
					varup = "g"+str(tr["down"])+str(tr["up"])+cur_metal+"_"+coll
					if(coll=="H2pa" or coll=="H2or"):
						MM[varM].append(varup+"*"+coll) #include colliders abundance if ortho or para H2
					else:
						MM[varM].append(varup+"*n(idx_"+coll+")") #include colliders abundance
					varMexist.append(varM)

				#build Mji as above
				for tr in transitions:
					varM = "M" + str(tr["up"]) + str(tr["down"]) + cur_metal
					real_variables.append(varM) #add to double variable list
					if(not(varM in varMexist)):
						try:
							MM[varM] = [Aijs[(tr["up"],tr["down"])]]
						except:
							MM[varM] = ["0e0"]
						MMij[varM] = [tr["up"],tr["down"]]
					coll = tr["coll"]
					varup = "g"+str(tr["up"])+str(tr["down"])+cur_metal+"_"+coll
					if(coll=="H2pa" or coll=="H2or"):
						MM[varM].append(varup+"*"+coll) #include colliders abundance if ortho or para H2
					else:
						MM[varM].append(varup+"*n(idx_"+coll+")") #include colliders abundance
					varMexist.append(varM)

				#write Mij and Mji definitions
				for k,v in MM.iteritems():
					full_cool += k + " = " + (" + &\n".join(v)) + "\n"
					full_cool += "\n"
				
				#prepares the matrix A and B
				full_cool += "!preparing matrix Ax=b\n"
				nlev = len(levels)
				Avar = "A"+cur_metal #A variable name
				real_variables.append(Avar+"("+str(nlev)+","+str(nlev)+")")
				Bvar = "B"+cur_metal #B variable name
				real_variables.append(Bvar+"("+str(nlev)+")")
				#build B and the first row of A which is 1
				full_cool += Bvar+"(:) = (/n(idx_"+cur_metal+"), " + (", ".join(["0d0"]*(nlev-1))) + "/)\n"
				full_cool += Avar+"(1,:) = (/" + (", ".join(["1d0"]*nlev)) + "/)\n"
				
				#build the other rows of A
				for ilev, lev in levels.iteritems():
					if(ilev==0): continue
					Arow = ["" for i in range(nlev)]
					for k,v in MMij.iteritems():
						if(ilev==v[0]): Arow[v[0]] += "-"+k
						if(ilev==v[1]): Arow[v[0]] += "+"+k
					Arow = ["0d0" if x=="" else x for x in Arow]
					full_cool += Avar+"("+str(ilev+1)+",:) = (/" + (", ".join(Arow)) + "/)\n"

				#include function to solve the linear system
				# depending on the number of levels
				full_cool += "\n"
				full_cool += "!solving Ax=b\n"
				if(nlev==2): 
					full_cool += "call mylin2("+Avar+"(:,:), "+Bvar+"(:))\n"
				elif(nlev==3): 
					full_cool += "call mylin3("+Avar+"(:,:), "+Bvar+"(:))\n"
				else:
					#LAPACK are called for more than 3 levels
					full_cool += "call mydgesv("+Avar+"(:,:), "+Bvar+"(:))\n"
					self.needLAPACK = True


				#prepares the cooling summing up the cooling from all the transitions
				full_cool += "\n"
				full_cool += "!computing cooling\n"
				cools = []
				trs = []
				for tr in transitions:
					if([tr["up"],tr["down"]] in trs): continue
					trs.append([tr["up"],tr["down"]])
					try:
						cool = Aijs[(tr["up"],tr["down"])] + "*"
						cool += str(float(levels[tr["up"]]["energy"]) - levels[tr["down"]]["energy"]) + " * "
						cool += Bvar + "(" + str(tr["up"]+1) + ")"
						cools.append(cool)
					except:
						pass
				full_cool += function_name + " = " + (" + &\n".join(cools)) + "\n"

				#insert the end of the function
				full_cool += "\n end function "+function_name+"\n"
	
				#replace the real variables needed by the cooling function
				vcool = ""
				uniq = []
				for x in real_variables:
					if(x in uniq): continue
					vcool += "real*8::"+x+"\n"
					uniq.append(x)
				#add ortho/para variables if needed
				if(needOrthoPara):
					vcool += "real*8::H2pa,H2or,o2p\n"

				full_cool = full_cool.replace("#KROME_replace_with_declarations", vcool)

				coolZ_functions.append([function_name, full_cool])

				continue
			
			#if block for rates, e.g. if(Tgas > 5d3)
			if(srow[:2]=="if"):
				ifcond, ifrate = srow.split(":") #if condition, and rate
				excitation_rates.append(ifcond.strip() +" "+ var_excitation +" = "+ifrate.strip()) #append rate
				krates.append(ifcond.strip() +" k("+ str(nkrates) +") = "+ifrate.strip())
				continue
			
			#read para/ortho
			if("para/ortho" in srow.replace(" ","")):
				aop = srow.split(":")[1].split("/")
				para2ortho = float(aop[0])/float(aop[1])
				ortho2para = float(aop[1])/float(aop[0])
				needOrthoPara = True
				print "Found para/ortho in "+fname+": "+str(para2ortho)

			#read ortho/para
			if("ortho/para" in srow.replace(" ","")):
				aop = srow.split(":")[1].split("/")
				para2ortho = float(aop[1])/float(aop[0])
				ortho2para = float(aop[0])/float(aop[1])
				needOrthoPara = True
				print "Found ortho/para in "+fname+": "+str(ortho2para)

			#read collider rates
			if(len(row.split(","))>=4):
				nkrates += 1
				#read collider, starting level, end level, rate F90 expression
				coll, lev_up, lev_down, rate = [x.strip() for x in srow.split(",",3)]
				coll = coll.replace("+","j").replace("-","k") #collider name for variable
				var_excitation = "g"+lev_up+lev_down+cur_metal+"_"+coll #excitation variable, g10C_Hp
				var_excitation_rev = "g"+lev_down+lev_up+cur_metal+"_"+coll #dexcitation variable, g01C_Hp

				if(not(coll in colliders)): colliders.append(coll) #append collider names
				real_variables.append(var_excitation) #add to double variable list
				excitation_rates.append(var_excitation +" = "+rate) #store excitation rates
				#print rate
				krates.append("k("+ str(nkrates) +") = "+rate+" !"+cur_metal+lev_up+lev_down+"_"+coll)
				kkrates.append(var_excitation+" = k("+str(nkrates) +")")
				transitions.append({"up":int(lev_up), "down":int(lev_down), "coll":coll}) #store level transition

			#read Aij
			if(len(row.split(","))==3):
				lup, ldown, Aij = [x.strip() for x in srow.split(",",4)]
				Aijs[(int(lup),int(ldown))] = Aij

		#copy local arrays and dictionaries to class attributes
		self.coolZ_functions = coolZ_functions
		self.coolZ_rates = krates
		self.coolZ_nkrates = nkrates
		self.coolZ_vars_cool = vars_cool
		
		#check for end statement								
		if(inmetal):
			print "ERROR: end metal statement missing in "+fname+"!"
			sys.exit()

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
		fh = open(self.srcFolder+"krome_commons.f90")
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
					fout.write("\tinteger,parameter::nrea=" + str(self.nrea) + "\n")
					fout.write("\tinteger,parameter::nmols=" + str(self.nmols) + "\n")
					fout.write("\tinteger,parameter::nspec=" + str(len(specs)) + "\n")
					fout.write("\tinteger,parameter::ndust=" + str(ndust) + "\n")
					fout.write("\tinteger,parameter::ndustTypes=" + str(self.dustTypesSize) + "\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
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
		constants.append(["n_mass","1.674920d-24","g"])
		constants.append(["ip_mass","1d0/p_mass","1/g"])
		constants.append(["clight","2.99792458e10","cm/s"]) 
		constants.append(["pi","3.14159265359d0","#"]) 
		constants.append(["eV_to_erg","1.60217646d-12","eV -> erg"]) 
		constants.append(["seconds_per_year","365d0*24d0*3600d0","yr -> s"]) 
		constants.append(["km_to_cm","1d5","km -> cm"]) 
		constants.append(["cm_to_Mpc","1.d0/3.08d24","cm -> Mpc"]) 
		constants.append(["kvgas_erg","8.d0*boltzmann_erg/pi/p_mass",""]) 
		constants.append(["pre_planck","2.d0*planck_erg/clight**2","erg/cm2*s3"]) 
		constants.append(["exp_planck","planck_erg / boltzmann_erg","s*K"]) 
		constants.append(["stefboltz_erg","5.670373d-5","erg/s/cm2/K4"])
		constants.append(["N_avogadro","6.0221d23","#"]) 
		constants.append(["Rgas_J","8.3144621d0","J/K/mol"]) 
		constants.append(["Rgas_kJ","8.3144621d-3","kJ/K/mol"])
		constants.append(["hubble","0.67d0","dimensionless"])
		constants.append(["Omega0","1.0d0","dimensionless"])
		constants.append(["Hubble0","1.d2*hubble*km_to_cm*cm_to_Mpc","1/s"])



		#********* CONSTANTS ****************
		fh = open(self.srcFolder+"krome_constants.f90")
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

			fh = open(self.srcFolder+"krome_user_commons.f90")
			fouta = open(self.buildFolder+"krome_user_commons.f90","w")

			for row in fh:
				row = row.replace("#KROME_header", get_licence_header(self.version, self.codename,self.shortHead))

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
		fh = open(self.srcFolder+"krome_subs.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_subs.f90","w")


		#create list of temperature shortcuts
		sclist = []
		for rea in reacts:
			sclist = get_Tshortcut(rea,sclist,coevars)

		#prepare shortcut definitions
		shortcutVars = ""
		for x in sclist:
			shortcutVars += "real*8::"+x.split("=")[0].strip()+"\n"

		#conserve
		krome_conserve = "" #init full string for the pragma replacement
		if(self.useConserve):
			skipa = ["CR","Tgas","dummy","g"] #skip these species
			multi = [] #species with shared atoms
			acount = dict() #store species per atom type (e.g. {"C":["C","CO","C2"], ...})
			has_multiple = False #check if spcies with shared atoms are present (e.g. CO but not H2)
			#loop on the species
			for x in specs:
				xname = x.name #species name
				if(xname in skipa): continue #cycle if the name is prensent in the skip list
				afound = 0 #cont founded type atoms
				#loop on the dictionary taht count the atoms in the species
				for a in x.atomcount2:
					if(a in ["+","-"]): continue #skip non-atoms
					afound +=1 #count founded atoms (for account of shared atoms)
					#append species to dictioanry if contains the atom a
					if(a in acount):
						acount[a].append(x)
					else:
						acount[a] = [x]
				if(afound>1): 
					has_multiple = True #flag for shared atoms
					multi.append(x.name)
			#if species with shared atoms warns the user (also in the subs file)
			if(has_multiple):
				krome_conserve += "!WARNING: found species with different atoms:\n"
				krome_conserve += "!conservation function may be non-accurate\n"
				krome_conserve += "!the approximation is valid when the following species\n"
				krome_conserve += "! are smaller compared to the others\n"
				krome_conserve += "! "+(", ".join(multi))+"\n"
				print
				print "WARNING: found species with different atoms:"
				print " conservation function may be non-accurate"
			print

			#loop on the found atoms. k=atom, v=list of species with k
			for (k,v) in acount.iteritems():
				if(k=="E"): continue #skip electrons
				aadd = [] #parts of the summation for ntot
				sdiff = "" #string for scaling
				#loop on the species
				for x in v:
					mult = (str(x.atomcount2[k])+"d0*" if x.atomcount2[k]>1 else "") #multiplication factor
					aadd.append(mult+"n("+x.fidx+")") #append species density with factor
					sdiff += "no("+x.fidx+") = n("+x.fidx+") * factor\n" #rescaling
				sadd = "ntot = " + (" &\n + ".join(aadd)) #current total density of the species k
				saddi = "nitot = " + (" &\n + ".join([y.replace("n(","ni(") for y in aadd])) #initial total density of the species k
				#prepare replacing string
				krome_conserve += "\n!********** "+k+" **********\n"
				krome_conserve += sadd + "\n"
				krome_conserve += saddi + "\n"
				krome_conserve += "factor = nitot/ntot\n"
				krome_conserve += sdiff + "\n"
				krome_conserve += "\n"
	
		nmax = 30
		if(len(specs)>nmax and self.useConserve):
			print "WARNING: more than "+str(nmax)+" species, -conserve disabled!"
			krome_conserve = "" #with more than NMAX species conserve only electrons

		has_electrons = False #check if electrons are present
		#check if electrons are present
		for x in specs:
			if(x.name=="E"): 
				has_electrons = True #check if electrons are present
				break

		#charge conservation
		if(has_electrons and self.useConserveE):
			consE = "n(idx_E) = max("
			for x in specs:
				if(x.name=="E"): continue #skip electron
				if(x.charge==0): continue #skip neutrals
				#prepare charge multiplicator
				mult = ""
				if(x.charge>0): mult = "+"
				if(x.charge>1): mult = "+" + str(x.charge) + "d0*"
				if(x.charge<0): mult = "-"
				if(x.charge<-1): mult = str(x.charge) + "d0*"
				consE += " &\n"+mult+"n("+x.fidx+")"
			consE += ", 1d-40)"
			krome_conserve += "\n!********** E **********\n"
			krome_conserve += consE + "\n"
		

		#loop on src file and replace pragmas
		skip = False
		for row in fh:

			srow = row.strip()

                        #skip when find IF pragmas
                        if(srow == "#IFKROME_useShieldingWG11" and not(self.useShieldingWG11)): skip = True
                        if(srow == "#IFKROME_useShieldingDB96" and not(self.useShieldingDB96)): skip = True
		        if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue #skip

			#replace the small value for rates according to the maximum number of products 
			if("#KROME_small" in srow):
				if(self.useTabs):
					fout.write(srow.replace("#KROME_small","0d0")+"\n")
					continue					
				maxprod = 0
				for x in reacts:
					maxprod = max(len(x.products),maxprod)
				mysmall = "1d-40/("+("*".join(["nmax"]*maxprod))+")"
				if(maxprod==0): mysmall = "0d0"
				fout.write(srow.replace("#KROME_small",mysmall)+"\n")
				continue

			#write reaction rates in coe function
			if(srow == "#KROME_krates"):
				for x in reacts:
					#build temperature limit IF
					sTlimit = ""
					hasTlim = (x.hasTlimitMin or x.hasTlimitMax) #Tmin or Tmax are present
					Tlimfound = False #flag to check if enfif is needed after the reaction rate
					if(x.kphrate==None and self.useTlimits and hasTlim):
						Tlimfound = True #need to close the if statement opened here
						sTlimit = "if("
						if(x.hasTlimitMin): sTlimit += "Tgas."+x.TminOp+"."+x.Tmin #Tmin is present
						if(x.hasTlimitMin and x.hasTlimitMax): sTlimit += " .and. " #Tmin and Tmax are present
						if(x.hasTlimitMax): sTlimit += "Tgas."+x.TmaxOp+"."+x.Tmax #Tmax is present
						sTlimit += ") then\n"
					kstr = "!" + x.verbatim+"\n" #reaction header
					kstr += "\t" + sTlimit + x.ifrate + " k("+str(x.idx)+") = " + x.krate #limit+extraif+rate
					if(Tlimfound): kstr += "\nend if" #close the if statement for temperature
					kstr = truncF90(kstr, 60,"*") #truncates long reaction rates
					fout.write(truncF90(kstr, 60,"/")+"\n\n") #truncate
			elif(srow == "#KROME_conserve"):
				fout.write(krome_conserve+"\n")
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
			elif(srow == "#KROME_imasses"):
				for x in specs:
					myimass = 0e0
					if(x.mass!=0e0): myimass = 1e0/x.mass
					imassrow = "\tget_imass("+str(x.idx)+") = " + str(myimass).replace("e","d") + "\t!" + x.name + "\n"
					fout.write(imassrow.replace("0.0","0.d0"))
			elif(srow == "#KROME_zatoms"):
				for x in specs:
					zatomrow = "\tget_zatoms("+str(x.idx)+") = " + str(x.zatom) + "\t!" + x.name + "\n"
					fout.write(zatomrow)
					
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
			elif(srow == "#KROME_qeff"):
				for x in reacts:
					sqeff = "\tget_qeff("+str(x.idx)+") = "+str(x.qeff)+" !" + x.verbatim
					fout.write(sqeff+"\n")
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
			elif(srow == "#KROME_shortcut_variables"):
				fout.write(shortcutVars)
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			elif(srow == "#KROME_gamma"):
				is_multiline = False #flag for multiline gamma
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

				elif(self.typeGamma=="EXACT" or self.typeGamma=="VIB" or self.typeGamma=="ROT"):
					#extends Omukai+Nishi1998 eqs.5,6,7
					# and Boley+2007 (+erratum!) eqs.2,3
					header = "real*8::Tgas,invTgas,x,expx,ysum,gsum,mosum,gvib\n"
					header += "real*8::Tgas_vib,invTgas_vib\n"
					gamma = "invTgas_vib = 1d0/Tgas_vib\n\n"
					#gamma += "nH = get_Hnuclei(n(:))\n\n"
					gi_vars = []
					g_vars = []
					di_vars = []
					mo_vars = []
					smallest_ve = 1e99
					print 
					for mol in specs:
						#monoatomic
						if(mol.natoms==1):
							gi_mo = "1.5d0" #where gi=1/(gamma_atom-1), where gamma_atom=5/3
							mo_vars.append(mol.fidx)
						#diatomic
						elif(mol.natoms==2):
							gtype = self.typeGamma
							#warning if vibrational constant not found
							if(mol.ve_vib=="__NONE__" and (gtype=="EXACT" or gtype=="VIB")):
								print "WARNING: no vibrational constant for "+mol.name+" in gamma calculation!"
							#warning if rotational constant not found
							if(mol.be_rot=="__NONE__" and (gtype=="EXACT" or gtype=="ROT")):
								print "WARNING: no rotational constant for "+mol.name+" in gamma calculation!"
							#continue if both constants were not found
							if(mol.ve_vib=="__NONE__" and mol.be_rot=="__NONE__"):
								continue
							gamma += "\n!evaluate 1/(gamma-1) for "+mol.name+"\n"
							#prepare the vibrational part
							vibpart = "0d0"
							if(mol.ve_vib!="__NONE__" and (gtype=="EXACT" or gtype=="VIB")):
								smallest_ve = min(smallest_ve, mol.ve_vib) #store the smallest vib constant
								xvar = "x = "+format_double(mol.ve_vib)+"*invTgas_vib\n"
								expvar = "expx = exp(x)\n"
								gamma += xvar
								gamma += expvar
								gamma += "gvib = 2d0*x*x*expx/(expx-1d0)**2\n"
								vibpart = "gvib"
							#prepare the rotational part
							rotpart = "2d0"
							if(mol.be_rot!="__NONE__" and (gtype=="EXACT" or gtype=="ROT")):
								if(mol.name=="H2"):
									rotpart = "gamma_rotop(Tgas, 3d0)"
								else:
									rotpart = "gamma_rot(Tgas, "+format_double(mol.be_rot)+")"
							di_vars.append(mol.fname)
							gi_vars.append("gi_"+mol.fname)
							gi = "gi_"+mol.fname+" = 0.5d0*(3d0 + "+rotpart+" + "+vibpart+")\n"
							gamma += gi

						#polyatomic
						else:
							pass #skip polyatomic
					#prepone variables declaration
					header += "real*8::"+(",".join(gi_vars)) + "\n"
					#write sums
					gamma += "\n!sum monotomic abundances\n"
					gamma += "mosum = " + (" + &\n".join(["n("+x+")" for x in mo_vars])) + "\n"
					gamma += "\n!sum all abundances\n"
					gamma += "ysum = mosum + "+(" + &\n".join(["n(idx_"+x+")" for x in di_vars]))+"\n"
					gamma += "\n!computes gamma\n"
					gamma += "gsum = mosum * "+gi_mo+" + "\
						+(" + &\n".join([("n(idx_"+x+")*gi_"+x) for x in di_vars]))+"\n"
					#add sum
					gamma += "krome_gamma = 1d0 + ysum/gsum\n"

					#append Tgas limit to avoid overflows on exp()
					header += "\n!avoid small Tgas that causes large x=a/Tgas below\n"
					header += "Tgas_vib = max(n(idx_Tgas), "+format_double(smallest_ve*1e-2) + ")\n"
					header += "Tgas = n(idx_Tgas)\n"
				
					#append gamma to the header
					gamma = header + gamma
					is_multiline = True
					
				else: 
					#user-defined gamma
					gamma = self.typeGamma

				#write gamma
				if(is_multiline):
					fout.write(gamma)
				else:
					fout.write("krome_gamma = " + gamma + "\n")
		
			else:
                                if(row[0]!="#"): fout.write(row)
		if(not(self.buildCompact)):
			fout.close()
		print "done!"


	##############################################
	def makePhoto(self):

		buildFolder = self.buildFolder
		reacts = self.reacts
		#********* PHOTO ****************
		print "- writing krome_photo.f90...",
		fh = open(self.srcFolder+"krome_photo.f90")
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
		coevars = self.coevars #copy coefficient varaibles
		#********* TABS ****************
		print "- writing krome_tabs.f90...",
		fh = open(self.srcFolder+"krome_tabs.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_tabs.f90","w")

		#include reactions that cannot be tabbed
		countNoTab = 0
		noTabReactions = ""
		sclist = [] #list of the temperature shortcuts
		for rea in self.reacts:
			if(not(rea.canUseTabs)): 
				noTabReactions += "coe_tab("+str(rea.idx)+") = "+rea.krate+"\n"
				countNoTab += 1
				sclist = get_Tshortcut(rea,sclist,coevars) #add shotcut if needed
			

		#if reactions that cannot be tabbed are found
		klist = kvars = ""
		if(countNoTab>0):
			if(len(coevars)!=0):
				#define variables
				kvars = "real*8::"+(",".join([x.strip() for x in coevars.keys()]))
				#variables initialization
				klist = [[k+" = "+v[1]+"\n",v[0]] for k,v in coevars.iteritems()] #this mess is to sort dict
				klist = sorted(klist, key=lambda x: x[1])
				klist = "".join([x[0] for x in klist])


		#replace pragmas
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
			row = row.replace("#KROME_logTup", "ktab_logTup = log10(min("+str(self.TmaxAuto)+",1d8))")
			row = row.replace("#KROME_define_vars",kvars)
			row = row.replace("#KROME_init_vars",klist)
			row = row.replace("#KROME_noTabReactions",noTabReactions)

			if(row.strip() == "#KROME_Tshortcuts"):
				ssc = ""
				for shortcut in sclist:
					ssc += shortcut + "\n"
				row = ssc

			if(self.useCustomCoe): row = row.replace("#KROMEREPLACE_customCoeFunction", self.customCoeFunction)

			if(len(row)==0): continue
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
		fh = open(self.srcFolder+"krome_dust.f90")
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
		fh = open(self.srcFolder+"krome_cooling.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_cooling.f90","w")


		#create coefficients, flux, and cooling functions for dH cooling
		i = 0
		dH_varsa = []
		dH_coe = dH_cool = ""
		if(self.useCoolingdH):
			idxs = []
			for rea in reacts:
				if(rea.idx in idxs): continue #skip reactions with the same index
				idxs.append(rea.idx)
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

			srow = row.strip()

			#cooling pragmas
			if(srow == "#IFKROME_useCoolingZ" and not(useCoolingZ)): skip = True
			if(srow == "#IFKROME_useCoolingdH" and (not(self.useCoolingdH) or len(dH_varsa)==0)): skip = True
			if(srow == "#IFKROME_useCoolingDust" and not(self.useCoolingDust)): skip = True
			if(srow == "#IFKROME_useCoolingAtomic" and not(self.useCoolingAtomic)): skip = True
			if(srow == "#IFKROME_useCoolingH2" and not(self.useCoolingH2)): skip = True
			if(srow == "#IFKROME_useCoolingH2GP" and not(self.useCoolingH2GP98)): skip = True
			if(srow == "#IFKROME_useCoolingHD" and not(self.useCoolingHD)): skip = True
			if(srow == "#IFKROME_useCoolingCompton" and not(self.useCoolingCompton)): skip = True
			if(srow == "#IFKROME_useCoolingExpansion" and not(self.useCoolingExpansion)): skip = True
			if(srow == "#IFKROME_useCoolingCIE" and not(self.useCoolingCIE)): skip = True
			if(srow == "#IFKROME_useCoolingContinuum" and not(self.useCoolingCont)): skip = True
			if(srow == "#IFKROME_useLAPACK" and not(self.needLAPACK)): skip = True #skip calls to LAPACK 

			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			#replace the small value for rates according to the maximum number of products 
			if("#KROME_small" in srow):
				if(self.useTabs):
					fout.write(srow.replace("#KROME_small","0d0")+"\n")
					continue					
				maxprod = 0
				for x in reacts:
					maxprod = max(len(x.products),maxprod)
				mysmall = "1d-40/("+("*".join(["nmax"]*maxprod))+")"
				if(maxprod==0): mysmall = "0d0"
				fout.write(srow.replace("#KROME_small",mysmall)+"\n")
				continue

			if(row.strip() == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			elif(row.strip() == "#KROME_nZrate"):
					fout.write("integer,parameter::nZrate="+str(self.coolZ_nkrates)+"\n")
			elif(row.strip() == "#KROME_coolingZ_call_functions"):
				for x in self.coolZ_functions:
					fout.write("cool = cool + "+x[0]+"(n(:),inTgas,k(:))\n")
			elif(row.strip() == "#KROME_coolingZ_custom_vars"):
				for x in self.coolZ_vars_cool:
					fout.write(x[0]+" = "+x[1]+"\n")
			elif(row.strip() == "#KROME_coolingZ_declare_custom_vars"):
				vcool = []
				for x in self.coolZ_vars_cool:
					vcool.append(x[0])
				if(len(vcool)>0): fout.write("real*8::"+(",".join(vcool))+"\n")
			elif(row.strip() == "#KROME_coolingZ_functions"):
				for x in self.coolZ_functions:
					fout.write(x[1]+"\n")
			elif(row.strip() == "#KROME_coolingZ_rates"):
				for x in self.coolZ_rates:
					fout.write(x+"\n")
			else:
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

		fh = open(self.srcFolder+"krome_heating.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_heating.f90","w")

		#create coefficients, flux, and cooling functions for dH cooling
		i = 0
		dH_varsa = []
		dH_coe = dH_heat = ""
		idxs = []
		if(self.useHeatingdH):
			for rea in reacts:
				if(rea.idx in idxs): continue #skip reactions with the same index
				idxs.append(rea.idx)
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
		if(self.useHeatingChem or self.useCoolingChem or self.useCoolingDISS):
			RPK = []
			#RPK is the list of the heating/cooling processes as 
			# [product_list, reactant_list, fortran_rate, heating/cooling_flag]
			if(self.useHeatingChem):
				RPK.append([["H","H","H"], ["H2","H"], "4.48d0*h2heatfac","H"])
				RPK.append([["H2","H","H"], ["H2","H2"], "4.48d0*h2heatfac","H"])
				RPK.append([["H-","H"], ["H2","E"], "3.53d0*h2heatfac","H"])
				RPK.append([["H2+","H"], ["H2","H+"], "1.83d0*h2heatfac","H"])
			if(self.useCoolingChem):
				RPK.append([["H","E"], ["H+","E","E"], "-13.6d0","C"])
				RPK.append([["HE","E"], ["HE+","E","E"], "-24.6d0","C"])
				RPK.append([["HE+","E"], ["HE++","E","E"], "-79.d0","C"])
			if(self.useCoolingChem or self.useCoolingDISS):
				RPK.append([["H2","H"], ["H","H","H"], "-4.48d0","C"])
				RPK.append([["H2","E"], ["H","H","E"], "-4.48d0","C"])
				RPK.append([["H2","H2"], ["H2","H","H"], "-4.48d0","C"])

			Rref = []
			Pref = []
			kref = []
			href = []
			for rpk in RPK:
				Rref.append(sorted(rpk[0])) #list of the reactants
				Pref.append(sorted(rpk[1])) #list fo the products
				kref.append(rpk[2]) #rate coefficient
				href.append(rpk[3]) #heating/cooling flag

			idxs = []
			for rea in reacts:
				if(rea.idx in idxs): continue #skip reactions with the same idx
				idxs.append(rea.idx)
				R = sorted([x.name for x in rea.reactants])
				P = sorted([x.name for x in rea.products])
				rmult = ("*".join(["n("+x.fidx+")" for x in rea.reactants]))
				for i in range(len(Rref)):
					if(Rref[i]==R and Pref[i]==P):
						headchem = "!"+rea.verbatim + " ("+("heating" if href[i]=="H"  else "cooling") + ")\n"
						hasTlim = (rea.hasTlimitMin or rea.hasTlimitMax) #Tmin or Tmax are present
						tklim = ""
						if(self.useTlimits and hasTlim):
							tklim = "if("
							if(rea.hasTlimitMin): tklim += "Tgas." + rea.TminOp + "."  + rea.Tmin #Tmin is present
							if(rea.hasTlimitMin and rea.hasTlimitMax): tklim += " .and. " #Tmax and Tmin are present
							if(rea.hasTlimitMax): tklim += "Tgas." + rea.TmaxOp + "." + rea.Tmax #Tmax is present
							tklim += ") then\n"
						HChem += headchem + tklim + "HChem = HChem + k("+str(rea.idx)+") * ("+kref[i] + "*"+rmult+")\n"
						if(self.useTlimits and hasTlim): HChem += "end if\n\n"
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

				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			else:
				if(row.strip() == "#IFKROME_useHeatingdH" and (not(self.useHeatingdH) or len(dH_varsa)==0)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingCompress" and not(self.useHeatingCompress)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingPhoto" and not(self.useHeatingPhoto)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				skipBool = (not(self.useHeatingChem) and not(self.useCoolingChem) and not(self.useCoolingDISS))
				if(row.strip() == "#IFKROME_useHeatingChem" and skipBool): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(skip): continue

				#replace the small value for rates according to the maximum number of products 
				if("#KROME_small" in row):
					if(self.useTabs):
						fout.write(row.replace("#KROME_small","0d0")+"\n")
						continue					
					maxprod = 0
					for x in reacts:
						maxprod = max(len(x.products),maxprod)
					mysmall = "1d-40/("+("*".join(["nmax"]*maxprod))+")"
					if(maxprod==0): mysmall = "0d0"
					row = row.replace("#KROME_small",mysmall)

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
		#dnw_grouped = self.dnw_grouped
		dnw = self.dnw
		neq = len(specs)
		solver_MF = self.solver_MF
		Tgas_species = self.Tgas_species
		coevarsODE = self.coevarsODE

		#*********ODE****************
		#write parameters in krome_ode.f90
		print "- writing krome_ode.f90...",

		fh = open(self.srcFolder+"krome_ode.f90")

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

						#add init flux var
						if(not(self.humanFlux)): 
							for x in self.reacts:
								fout.write("kflux("+str(x.idx)+") = "+x.RHS+"\n")
							fout.write("\n")
					
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
			elif(srow == "#KROME_flux_variables" and not(self.humanFlux)):
				#add var declaration for flux
				#for x in self.reacts:
				#	fout.write("real*8::"+x.RHSvar+"\n")
				fout.write("real*8::kflux("+str(len(self.reacts))+")\n")
				fout.write("\n")

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
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
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
			elif(srow == "#KROME_JAC_PDX"):
				if(not(self.doJacobian)): continue
				spdj = ""
				for i in range(neq):
					speci = specs[i]
					for j in range(neq):
						specj = specs[j]
						if(self.jsparse[j][i]==1):
							org = "pdj("+str(j+1)+")"
							rep = "pd("+str(j+1)+","+str(i+1)+")"
							orgT = "pdj(idx_Tgas)"
							repT = "pd(idx_Tgas,"+str(i+1)+")"
							spdj += "!d["+str(specj.name)+"_dot]/d["+str(speci.name)+"]\n"
							spdj += self.jac[j][i].replace(org,rep).replace(orgT,repT)+"\n\n"
							
				fout.write(spdj)

			elif(srow == "#KROME_JAC_PD"):
				if(not(self.doJacobian)): continue
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
		fh = open(self.srcFolder+"krome_user.f90")
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

			if(srow == "#IFKROME_useStars" and not(self.useStars)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_use_cooling" and not(self.use_cooling)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_use_thermo" and not(self.use_thermo)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			if(srow == "#KROME_species"):
				for x in specs:
					fout.write("\tinteger,parameter::" + "KROME_"+x.fidx + " = " + str(x.idx) +"\t!"+x.name+"\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			elif(srow == "#KROME_zero_electrons"):
				#check if electron exists
				for x in specs:
					if(x.name=="E"):
						fout.write("x(idx_e) = 0.d0\n")
						break
			elif(srow == "#KROME_electrons_balance"):
				fout.write("x(idx_e) = ee\n")
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
				fout.write("\tinteger,parameter::krome_nrea=" + str(self.nrea) + "\n")
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
		fh = open(self.srcFolder+"krome_reduction.f90")

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
				#print self.maxnreag
				if(self.maxnreag>0):
					fout.write("arr_flux(i) = k(i)*"+("*".join(["n(r"+str(j+1)+")" for j in range(self.maxnreag)]))+"\n")

			if(row[0]!="#"): fout.write(row)
		if(not(self.buildCompact)):
			fout.close()

		print "done!"

	##############################
	def makeStars(self):		
		buildFolder = self.buildFolder
		#********* STARS ****************
		#intended for nuclear networks of stars
		print "- writing krome_stars.f90...",
		fh = open(self.srcFolder+"krome_stars.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_stars.f90","w")


		skip = False
		for row in fh:
			srow = row.strip()

			if(srow == "#IFKROME_useStars" and not(self.useStars)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue
			#computes screening for 3body reactions (as screen products [1,2]*[1+2,3])
			if(srow == "#KROME_stars_3body"):
				stars3Body = "" #3body string to replace pragma
				idxs = []
				for rea in self.reacts:
					if(rea.idx in idxs): continue
					idxs.append(rea.idx)
					if(len(rea.reactants)==3 and not(True in rea.curlyR)):
						sdx = str(rea.idx) #string index
						stars3Body += "\n!3body: "+rea.verbatim+"\n"
						stars3Body += "z12 = zz(arr_r1("+sdx+")) + zz(arr_r2("+sdx+"))\n"
						stars3Body += "scr12 = stars_screen(Tgas,rho,n(:), zz(arr_r1("+sdx+")), zz(arr_r2("+sdx+")))\n"
						stars3Body += "scr23 = stars_screen(Tgas,rho,n(:), z12, zz(arr_r3("+sdx+")))\n"
						stars3Body += "k("+sdx+") = ko("+sdx+") * scr12 * scr23 \n"
				fout.write(stars3Body)
			elif(srow=="#KROME_stars_energy"):
				starsE = "" #energy string replacing pragma
				idxs = []
				for rea in self.reacts:
					if(rea.idx in idxs): continue #avoid reactions with same index
					idxs.append(rea.idx)
					sdx = str(rea.idx) #string index
					starsE += "flux("+sdx+") = "+rea.RHS+" !"+rea.verbatim+"\n"
				fout.write(starsE)

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
			fh = open(self.srcFolder+"kromeF90.f90")
		else:
			fh = open(self.srcFolder+"krome.f90")

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

			if(srow == "#IFKROME_conserve" and not(self.useConserve) and not(self.useConserveE)): skip = True
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

			if(srow == "#IFKROME_useStars" and not(self.useStars)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_useCoolingZ" and not(self.useCoolingZ)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_ierr" and not(self.useIERR)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			if(srow == "#IFKROME_noierr" and (self.useIERR)): skip = True
			if(srow == "#ENDIFKROME"): skip = False

			ierr = ""
			if(self.useIERR): ierr = ",ierr"
			if(self.useDust):
				row = row.replace("#KROME_dust_arguments",",xdust"+ierr)
			else:
				row = row.replace("#KROME_dust_arguments",""+ierr)

			row = row.replace("#KROME_ATOL",str(ATOL))
			row = row.replace("#KROME_RTOL",str(RTOL))

			if(skip): continue

			if(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
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
			elif(srow == "#KROME_iaja_parameters"):
				fout.write("integer,parameter::niauser="+str(len(self.ia))+",njauser="+str(len(self.ja))+"\n")
			elif(srow == "#KROME_maxord" and self.maxord!=0):
				fout.write("iopt = 1 !activate optional inputs\n")
				fout.write("IWORK(5) = "+str(self.maxord)+" !maximum integration order\n")
			elif(srow == "#KROME_MF"):
				fout.write("MF = "+str(self.solver_MF)+"\n")
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
		print "- copying optical data for dust..."
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
		print "- copying solver to /build..."
		if(self.useDvodeF90):
			shutil.copyfile("solver/dvode_f90_m.f90", buildFolder+"dvode_f90_m.f90")
		else:
			shutil.copyfile("solver/opkdmain.f", buildFolder+"opkdmain.f")
			shutil.copyfile("solver/opkda1.f", buildFolder+"opkda1.f")
			shutil.copyfile("solver/opkda2.f", buildFolder+"opkda2.f")


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
				indentF90(buildFolder+"krome_stars.f90")
				indentF90(buildFolder+"krome_reduction.f90")
				indentF90(buildFolder+"krome_subs.f90")
				indentF90(buildFolder+"krome_tabs.f90")
				indentF90(buildFolder+"krome_user.f90")

		print "done!"
	#########################################
	#copy fsrc to fout file and replace the list in pragmas with the list in repls.
	# if trim the each line is trimmed and a return \n is added when write it again
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
		if(not(os.path.exists(ramsesFolder))): os.makedirs(ramsesFolder)
		specs = self.specs
		ramses_offset = str(self.ramses_offset)

		#some initial abundances. if not found set default
		# all in 1/cm3, except for T which is K
		ndef = {"H": 7.5615e-1,
			"E": 4.4983e-8,
			"H+": 8.1967e-5,
			"HE": 2.4375e-1,
			"H2": 1.5123e-6,
			"Tgas" : 200,
			"default":1e-40
		}

		excl = ["CR","g","Tgas","dummy"] #avoid specials

		#count species excluding what is conteinted in excl list
		chemCount = 0
		for x in specs:
			if(x.name in excl): continue
			chemCount += 1

		#amr_parameters
		#just copy the file fname to the build/ramses folder
		fname = "amr_parameters.f90"
		self.replacein(pfold+fname, ramsesFolder+fname, ["aaa"], ["aaa"])
		indentF90(ramsesFolder+fname)

		#condinit
		#prepares the initial conditions and copy fname
		cheminit = " q(1:nn,ndim+3) = "+str(ndef["Tgas"])+"     !Set temperature in K\n"
		ichem = 0
		fname = "condinit.f90"
		#loop on species
		for x in specs:
			#skip species in exl list
			if(x.name in excl): continue
			ichem += 1
			#check if species is in init array (ndef) else default
			if(x.name in ndef):
				sdef = str(ndef[x.name]) #default value from array
			else:
				sdef = str(ndef["default"]) #default values if not present in array
			cheminit += "q(1:nn,ndim+3+"+str(ichem)+")  = "+sdef+"  !"+x.name+"\n"
		#replace initialization
		self.replacein(pfold+fname, ramsesFolder+fname, ["#KROME_init_chem"], [cheminit])
		indentF90(ramsesFolder+fname)

		#cooling_fine
		#prepare the array for krome and back (ramses->krome->ramses)
		# updateueq: copy from ramses 2dim array to 1dim array for krome (unoneq<-uold)
		# scaleueq: scale array from code units (RAMSES) to 1/cm3 (KROME) (unoneq<-unoneq)
		# bkscaleueq: scale array from 1/cm3 (KROME) to code units (RAMSES) (unoneq->unoneq)
		# bkupdateueq: copy from 1dim array of krome to 2dim array of ramses (unoneq->uold)
		updateueq = scaleueq = bkscaleueq = bkupdateueq = ""
		ichem = 0
		fname = "cooling_fine.f90"
		for x in specs:
			ichem += 1
			if(not(x.name in excl)):
				updateueq += "unoneq("+str(ichem)+") = uold(ind_leaf(i),ndim+"+ramses_offset+"+"+str(ichem)+") !"+x.name+"\n"
				if(x.mass>0e0): scaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*scale_d/"+str(x.mass)+" !"+x.name+"\n"
				bkscaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*"+str(x.mass)+"/scale_d !"+x.name+"\n"
				bkupdateueq += "uold(ind_leaf(i),ndim+"+ramses_offset+"+"+str(ichem)+") = unoneq("+str(ichem)+")\n"
		org = ["#KROME_update_unoneq","#KROME_scale_unoneq","#KROME_backscale_unoneq","#KROME_backupdate_unoneq"]
		new = [updateueq, scaleueq, bkscaleueq, bkupdateueq]
		#replace pragmas (org) with expressions (new)
		self.replacein(pfold+fname, ramsesFolder+fname, org, new)
		indentF90(ramsesFolder+fname)

		#cooling_module
		# simply copy the cooling_module into build/ramses
		fname = "cooling_module.f90"
		self.replacein(pfold+fname,ramsesFolder+fname, [], [])
		indentF90(ramsesFolder+fname)

		#hydro_parameters
		# simply copy the hydro_parameters into build/ramses
		# extend nvar according to KROME species
		fname = "hydro_parameters.f90"
		self.replacein(pfold+fname, ramsesFolder+fname, [], [])
		#indentF90(ramsesFolder+fname)

		#init_flow_fine
		fname = "init_flow_fine.f90"
		init_array = "if(ivar==ndim+"+ramses_offset+")  init_array = 1.356d-2/aexp**2 ! T in K\n"
		ichem = 0
		for x in specs:
			if(x.name in excl): continue
			ichem += 1
			#check if species is contained in the ndef array (see above)
			if(x.name in ndef):
				sdef = str(ndef[x.name]) #default value from array
			else:
				sdef = str(ndef["default"]) #default value if not present in array
			init_array += "if(ivar==ndim+"+ramses_offset+"+"+str(ichem)+")  init_array = "+sdef+"  !"+x.name+"\n"
		#replace pragma and copy the file to the build/ramses
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
		self.replacein(pfold+fname,buildFolder+fname,["#KROME_nvar"],["#this must be NDIM+"+str(chemCount+3)], False)

		#move the krome files into the ramses patch folder
		shutil.move(buildFolder+"krome_all.f90", ramsesFolder+"krome_all.f90")
		shutil.move(buildFolder+"krome_user_commons.f90", ramsesFolder+"krome_user_commons.f90")
		shutil.move(buildFolder+"opkda1.f", ramsesFolder+"opkda1.f")
		shutil.move(buildFolder+"opkda2.f", ramsesFolder+"opkda2.f")
		shutil.move(buildFolder+"opkdmain.f", ramsesFolder+"opkdmain.f")


	#########################################
	def ramsesTH_patch(self):
		pfold = "patches/ramsesTH/" #source folder
		ramsesFolder = self.buildFolder+"krome_ramsesTH_patch/" #destination folder
		buildFolder = self.buildFolder
		if(not(os.path.exists(ramsesFolder))): os.makedirs(ramsesFolder)
		specs = self.specs

		#some initial abundances. if not found set default
		# all in 1/cm3, except for T which is K
		ndef = {"H": 7.5615e-1,
			"E": 4.4983e-8,
			"H+": 8.1967e-5,
			"HE": 2.4375e-1,
			"H2": 1.5123e-6,
			"O": 0.003792e0,
			"C": 0.001269e0,
			"Tgas" : 200e0,
			"default":1e-40
		}

		excl = ["CR","g","Tgas","dummy"] #avoid specials

		#count species excluding what is conteinted in excl list
		chemCount = 0
		for x in specs:
			if(x.name in excl): continue
			chemCount += 1


		#cooling_fine
		#prepare the array for krome and back (ramses->krome->ramses)
		# updateueq: copy from ramses 2dim array to 1dim array for krome (unoneq<-uold)
		# scaleueq: scale array from code units (RAMSES) to 1/cm3 (KROME) (unoneq<-unoneq)
		# bkscaleueq: scale array from 1/cm3 (KROME) to code units (RAMSES) (unoneq->unoneq)
		# bkupdateueq: copy from 1dim array of krome to 2dim array of ramses (unoneq->uold)
		updateueq = scaleueq = bkscaleueq = bkupdateueq = ""
		ichem = 0
		fname = "cooling_fine.f90"
		for x in specs:
			ichem += 1
			if(not(x.name in excl)):
				updateueq += "unoneq("+str(ichem)+") = uold(ind_leaf(i),ichem+1+"+str(ichem-1)+") !"+x.name+"\n"
				if(x.mass>0e0): scaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*scale_d/"+str(x.mass)+" !"+x.name+"\n"
				bkscaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*"+str(x.mass)+"*iscale_d !"+x.name+"\n"
				bkupdateueq += "uold(ind_leaf(i),ichem+1+"+str(ichem-1)+") = unoneq("+str(ichem)+") !"+x.name+"\n"
		org = ["#KROME_update_unoneq","#KROME_scale_unoneq","#KROME_backscale_unoneq","#KROME_backupdate_unoneq"]
		new = [updateueq, scaleueq, bkscaleueq, bkupdateueq]
		#replace pragmas (org) with expressions (new)
		self.replacein(pfold+fname, ramsesFolder+fname, org, new)
		indentF90(ramsesFolder+fname)


		#prepares abundances.nml
		abnml = "!This file contains the initialization for the species\n"
		abnml += "! employed by KROME. Change them according to your needs\n"
		ichem = 0
		comma = "" #add a separator here for namelist
		for x in specs:
			if(chemCount==ichem+1): comma = ""
			if(x.name in excl): continue
			#check if species is contained in the ndef array (see above)
			if(x.name in ndef):
				sdef = format_double(ndef[x.name]) #default value from array
			else:
				sdef = format_double(ndef["default"]) #default value if not present in array
			if(ichem==0):
				abpart = "metal_region(1,1) = "
				absize = len(abpart)
				abpart += sdef+comma+(" "*(20-len(sdef)-len(comma)))
			else:
				abpart = (" "*absize)+sdef+comma+(" "*(20-len(sdef)-len(comma)))
			abnml += abpart+"!"+str(ichem+1)+": "+x.name+"\n"
			ichem += 1

		#write abundances.nml
		fh = open(buildFolder+"abundances.nml","w")
		fh.write(abnml)
		fh.close()

		#copy cooling
		fname = "cooling.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		indentF90(ramsesFolder+fname)

		#copy Makefile
		fname = "Makefile"
		liblapack = ""
		if(self.needLAPACK): liblapack = "LIBS += -llapack"
		#+1 after ichem is for adiabatic index
		self.replacein(pfold+fname, ramsesFolder+fname, ["#KROME_NMOLS","#KROME_useLAPACK"], [str(ichem+1),liblapack], False)

		#copy Makefile.dep
		fname = "Makefile.dep"
		shutil.copy(pfold+fname,ramsesFolder+fname)


		#move the krome files into the ramses patch folder
		shutil.move(buildFolder+"krome_all.f90", ramsesFolder+"krome_all.f90")
		shutil.move(buildFolder+"krome_user_commons.f90", ramsesFolder+"krome_user_commons.f90")
		shutil.move(buildFolder+"opkda1.f", ramsesFolder+"opkda1.f")
		shutil.move(buildFolder+"opkda2.f", ramsesFolder+"opkda2.f")
		shutil.move(buildFolder+"opkdmain.f", ramsesFolder+"opkdmain.f")



	###########################################
	def flash_patch(self):
		specs = self.specs


		#some initial fractions. if not found set default
		ndef = {"H": 0.76e0,
			"HE": 0.24e0,
			"H+": 1.2e-5,
			"H-": 2e-9,
			"He+": 1e-14,
			"He++": 1e-17,
			"H2": 2e-20,
			"default":1e-40
		}

		buildFolder = self.buildFolder
		flashFolder = buildFolder+"krome_flash_patch/"
		patchFolder = "patches/flash/"

		#create folders
		dirs = [flashFolder + "Driver/DriverMain"]
		dirs.append(flashFolder + "physics/sourceTerms/KromeChemistry/KromeChemistryMain")
		dirs.append(flashFolder + "Simulation/SimulationComposition/KromeChemistry")
		dirs.append(flashFolder + "Simulation/SimulationMain/Chemistry_Krome_Collapse")
		for pdir in dirs:
			if(not(os.path.exists(pdir))): os.makedirs(pdir)

		excl = ["CR","g","Tgas","dummy"] #species to exclude

		########Driver#############
		pFolder = "Driver/DriverMain/"
		flist = ["Driver_finalizeSourceTerms.F90","Driver_initSourceTerms.F90","Driver_sourceTerms.F90"]
		for fle in flist:
			shutil.copy(patchFolder+pFolder+fle, flashFolder+pFolder+fle)


		########physics#############
		#*******Config**********
		pFolder = "physics/sourceTerms/KromeChemistry/KromeChemistryMain/"
		species = ""
		speciesCount = 0
		for x in specs:
			if(x.name in excl): continue
			name = x.name.upper().replace("+","P").replace("-","M")
			if(name=="E"): name="ELEC"
			speciesCount += 1
			species += "SPECIES "+name+"\n"

		fname = "Config"
		self.replacein(patchFolder+pFolder+fname, flashFolder+pFolder+fname,["#KROME_spec_data"],[species])

		#*******build->physics****************
		pFolder = "physics/sourceTerms/KromeChemistry/KromeChemistryMain/"
		kromeFileList = ["krome_all.f90", "krome_user_commons.f90", "opkda1.f", "opkda2.f", "opkdmain.f"]
		for fl in kromeFileList:
			shutil.move(buildFolder+fl, flashFolder+pFolder+fl)

		#***********physics->physics***********
		pFolder = "physics/sourceTerms/KromeChemistry/"
		kromeFileList = ["KromeChemistry.F90","KromeChemistry_finalize.F90", 
			"KromeChemistry_init.F90", "KromeChemistry_interface.F90", "Makefile"]
		for fl in kromeFileList:
			shutil.copy(patchFolder+pFolder+fl, flashFolder+pFolder+fl)

		#***********physics_main->physics_main***********
		kromeFileList = ["KromeChemistry.F90", "KromeChemistry_data.F90", "KromeChemistry_init.F90", "Makefile"]
		pFolder = "physics/sourceTerms/KromeChemistry/KromeChemistryMain/"
		for fl in kromeFileList:
			shutil.copy(patchFolder+pFolder+fl, flashFolder+pFolder+fl)

		#**************pchem_mapNetworkToSpecies*********
		pFolder = "physics/sourceTerms/KromeChemistry/KromeChemistryMain/"
		fname = "pchem_mapNetworkToSpecies.F90"
		species = ""
		for x in specs:
			if(x.name in excl): continue
			name = x.name.upper().replace("+","P").replace("-","M")
			if(name=="E"): name="ELEC"
			species += "\tcase(\""+(x.name)+"\")\n"
			species += "\t\tspecieOut = "+name+"_SPEC\n"

		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_cases"],[species])
		indentF90(flashFolder+pFolder+fname)

		########Simulation#############
		#***********pfold->treeFolder***********
		kromeFileList = ["Config"]
		pFolder = "Simulation/SimulationComposition/KromeChemistry/"
		for fl in kromeFileList:
			shutil.copy(patchFolder+pFolder+fl, flashFolder+pFolder+fl)

		#******simulation_initSpecies
		fname = "Simulation_initSpecies.F90"
		pFolder = "Simulation/SimulationComposition/KromeChemistry/"
		self.replacein(patchFolder+pFolder+fname, flashFolder+pFolder+fname,["#KROME_specnum"],[str(speciesCount)])
		indentF90(flashFolder+pFolder+fname)

		#SpeciesList.txt
		fname = "SpeciesList.txt"
		pFolder = "Simulation/SimulationComposition/KromeChemistry/"
		spec_data = ""
		all_parts = []
		for x in specs:
			if(x.name in excl): continue
			name = x.name.upper().replace("+","P").replace("-","M")
			if(name=="E"): name="ELEC"
			#prepares gamma depending on the model employed for gamma
			if(self.typeGamma=="FULL"):
				if(x.is_atom): #monoatomic
					gamma = 5./3.
				else: #molecule
					gamma = 7./5.
			elif(self.typeGamma=="DEFAULT"):
				gamma = "1.66666666667d0"
			else:
				gamma = self.typeGamma
			all_parts.append([name, x.zatom, x.mass, x.neutrons, x.zatom-x.charge,gamma])
		all_parts = sorted(all_parts,key=lambda x:x[1]) #sort by atomic number
		for parts in all_parts:
			spec_data += ("".join([str(y)+(20-len(str(y)))*" " for y in parts]))	+ "\n"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_spec_data"],[spec_data])

		#************####Collapse example###**************
		pFolder = "Simulation/SimulationMain/Chemistry_Krome_Collapse/"
		specs_config = specs_par = specs_data =  specs_init = ""
		specs_block_vars = specs_block_if = specs_block_prop = ""
		ablocke = []
		for x in specs:
			if(x.name in excl): continue
			name = x.name.upper().replace("+","P").replace("-","M")
			if(name=="E"): name="ELEC"
			if(x.name in ndef):
				nn = ndef[x.name]
			else:
				nn = ndef["default"]
			parts = ["PARAMETER", "sim_x"+name, "REAL", nn]
			specs_config += ("".join([str(y)+(20-len(str(y)))*" " for y in parts]))	+ "\n"
			parts = ["sim_x"+name, "=", nn]
			specs_par += ("".join([str(y)+(20-len(str(y)))*" " for y in parts]))	+ "\n"
			specs_data += "real, save :: sim_x"+name +"\n"
			specs_init += "call RuntimeParameters_get(\"sim_x"+name+"\", sim_x"+name+")\n"
			specs_block_vars += "real :: "+name.lower()+"A\n"
			specs_block_if += "if("+name+"_SPEC > 0) massFraction("+name+"_SPEC) = max(sim_x"+name+", smallx)\n"
			specs_block_prop += "call Multispecies_getProperty("+name+"_SPEC,A,"+name.lower()+"A)\n"
			if(x.charge!=0 and x.name!="E"):
				if(x.charge==1):
					emult = " +"
				elif(x.charge==-1):
					emult = " -"
				else:
					emult = " "+str(int(x.charge))+"*"
					if(x.charge>0): emult = "+" + emult
				ablocke.append(emult+"massFraction("+name+"_SPEC)/"+name.lower()+"A")
		specs_block_e = "if (ELEC_SPEC > 0) massFraction(ELEC_SPEC) = max( elecA*(&\n"
		specs_block_e += ("&\n".join(ablocke))
		specs_block_e += "&\n), smallx)"

		fname = "Config"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_specs_config"],[specs_config])
		indentF90(flashFolder+pFolder+fname)

		fname = "flash.par"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_specs_par"],[specs_par])
		indentF90(flashFolder+pFolder+fname)

		fname = "flash.par_1d"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_specs_par"],[specs_par])
		indentF90(flashFolder+pFolder+fname)

		fname = "Simulation_data.F90"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_specs_data"],[specs_data])
		indentF90(flashFolder+pFolder+fname)

		fname = "Simulation_init.F90"
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,["#KROME_specs_init"],[specs_init])
		indentF90(flashFolder+pFolder+fname)

		fname = "Simulation_initBlock.F90"
		pragmas = ["#KROME_specs_block_vars", "#KROME_specs_block_if","#KROME_specs_block_prop","#KROME_specs_block_e"]
		reps = [specs_block_vars, specs_block_if,specs_block_prop,specs_block_e]
		self.replacein(patchFolder+pFolder+fname,flashFolder+pFolder+fname,pragmas,reps)
		indentF90(flashFolder+pFolder+fname)

		#***********pfold->treeFolder***********
		kromeFileList = ["Makefile"]
		pFolder = "Simulation/SimulationMain/Chemistry_Krome_Collapse/"
		for fl in kromeFileList:
			shutil.copy(patchFolder+pFolder+fl, flashFolder+pFolder+fl)

	########################################
	#cut the arg sting using sep as separtor
	# when a piece is longer than largmax append the
	# chracter rep
	def linebreakerC(self,arg,sep,largmax = 40,rep="\n"):
		
		aarg = arg.split(sep)
		sarg = ""
		larg = 0
		for x in aarg[:len(aarg)-1]:
			sarg += x + sep
			larg += len(x+sep)
			if(larg>largmax):
				sarg += rep
				larg = 0
		sarg += aarg[len(aarg)-1]
		return sarg

	###########################################
	def enzo_patch(self):
		buildFolder = self.buildFolder
		enzoFolder = buildFolder+"krome_enzo_patch/"
		patchFolder = "patches/enzo/"

		if(not(os.path.exists(enzoFolder))): os.makedirs(enzoFolder)

		specs = self.specs
		excl = ["CR","g","Tgas","dummy"] #species to exclude
		speciesCount = 0
		krome_driver_args = krome_driver_rprec = krome_driver_scale = ""
		krome_driver_minval = krome_driver_dom = krome_driver_mod = ""
		krome_identify_num = ""
		krome_solve_args = krome_solve_baryon = ""
		krome_driver_suma = []
		kdrive_args = []
		krome_identify_identifya = []
		krome_identify_zeroa = []
		krome_identify_binarya = []
		krome_identify_vfail1a = []
		krome_identify_vfail2a = []
		krome_solve_numa = []
		krome_solve_identifya = []
		krome_grid_identifya = []
		for x in specs:
			if(x.name in excl): continue
			uname = x.name.upper() #upper-case name
			uname = uname.replace("HE","He") #Helium is lowecase
			#enzo-like names, e.g. H- -> HM, and H+ -> HII
			if("-" in uname):
				name = uname.replace("-","M") #anions
			else:
				name = (uname+"I").replace("+","I") #neutral and ions
			extname = name+"Density"
			if(name=="EI"): 
				name = "De" #electron is special
				extname = "ElectronDensity"
			speciesCount += 1 #increases species count
			#1. DRIVER file pragama
			krome_driver_args += name+", " #function arguments
			if(speciesCount%4==0): krome_driver_args += "&\n"
			krome_driver_rprec += "real*8::"+name+"(in,jn,kn)\n"
			krome_driver_scale += name+"(i,j,k) = "+name+"(i,j,k) * factor\n"
			krome_driver_minval += name+"(i,j,k) = max("+name+"(i,j,k), krome_tiny)\n"
			dmult = ""
			xzn = x.zatom+x.neutrons
			if(xzn>1): dmult = "* "+str(1e0/xzn)+"d0"
			krome_driver_dom += "krome_x(krome_"+x.fidx+") = "+name+"(i,j,k) * dom "+dmult+"\n"
			krome_driver_suma.append("+"+name+"(i,j,k) "+dmult)
			dmult2 = ""
			if(xzn>1): dmult2 = "* "+str(xzn)+"d0"
			krome_driver_mod += name+"(i,j,k) = krome_x(krome_"+x.fidx+") * idom "+dmult2+"\n"

			#2. Grid_IdentifySpeciesFieldsKrome.C
			krome_identify_identifya.append("int &"+name+"Num")
			krome_identify_zeroa.append(name+"Num")
			krome_identify_binarya.append(name+"Num<0")
			krome_identify_vfail1a.append(name+"=%\"ISYM\"")
			krome_identify_vfail2a.append(name+"Num")
			krome_identify_num += " "+name+"Num = FindField("+extname+", FieldType, NumberOfBaryonFields);\n"

			#3. Grid_SolveRateAndCoolEquations
			krome_solve_args += "float *"+name+", "
			krome_solve_baryon += "BaryonField["+name+"Num], "
			krome_solve_numa.append(name+"Num")
			krome_solve_identifya.append(name+"Num")
			
			#4. Grid.h
			krome_grid_identifya.append("int &"+name+"Num")

		if(speciesCount%4!=0): krome_driver_args += "&"
		krome_driver_sum = (" &\n".join(krome_driver_suma))

		krome_identify_identify = self.linebreakerC((", ".join(krome_identify_identifya)), ",")
		krome_identify_zero = self.linebreakerC((" = ".join(krome_identify_zeroa)), "=")+" = 0;"
		krome_identify_binary = self.linebreakerC((" || ".join(krome_identify_binarya)), "||")
		krome_identify_vfail1 = "\""+(", ".join(krome_identify_vfail1a))+"\\n\""
		krome_identify_vfail2 = self.linebreakerC((", ".join(krome_identify_vfail2a)), ",")

		krome_solve_args = self.linebreakerC(krome_solve_args, ",")
		krome_solve_num = "int "+ self.linebreakerC((", ".join(krome_solve_numa)), ",")+";"
		krome_solve_identify = self.linebreakerC((", ".join(krome_solve_identifya)), ",")
		krome_solve_baryon = self.linebreakerC(krome_solve_baryon, ",")

		krome_grid_identify = self.linebreakerC((", ".join(krome_grid_identifya)), ",")

		#1. replace
		fname = "krome_driver.F90"
		prags = ["#KROME_args","#KROME_rprec","#KROME_scale","#KROME_minval","#KROME_dom","#KROME_mod","#KROME_sum"]
		reps = [krome_driver_args,krome_driver_rprec,krome_driver_scale,krome_driver_minval]
		reps += [krome_driver_dom,krome_driver_mod,krome_driver_sum]
		self.replacein(patchFolder+fname, enzoFolder+fname, prags, reps)
		indentF90(enzoFolder+fname)

		#2. replace in Grid_IdentifySpeciesFieldsKrome.C
		fname = "Grid_IdentifySpeciesFieldsKrome.C"
		prags = ["#KROME_identify", "#KROME_zero","#KROME_binary", "#KROME_vfail1", "#KROME_vfail2", "#KROME_num"]
		reps = [krome_identify_identify, krome_identify_zero, krome_identify_binary]
		reps += [krome_identify_vfail1, krome_identify_vfail2, krome_identify_num]
		self.replacein(patchFolder+fname, enzoFolder+fname, prags, reps, False)

		#3. replace in Grid_SolveRateAndCoolEquations.C
		fname = "Grid_SolveRateAndCoolEquations.C"
		prags = ["#KROME_args", "#KROME_num", "#KROME_identify","#KROME_baryon"]
		reps = [krome_solve_args, krome_solve_num, krome_solve_identify, krome_solve_baryon]
		self.replacein(patchFolder+fname, enzoFolder+fname, prags, reps, False)

		#4. replace in Grid.h
		fname = "Grid.h"
		self.replacein(patchFolder+fname, enzoFolder+fname, ["#KROME_identify"], [krome_grid_identify], False)

		#5. copy others
		fname = "evaluate_temp.F90"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "InitializeRateData.C"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "krome_initab.F90"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "Make.config.objects"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "kromebuild.sh"
		shutil.copy(patchFolder+fname, enzoFolder+fname)

		#6. move others
		flist = ["krome_all", "krome_user_commons"]
		for fle in flist:
			shutil.move(buildFolder+fle+".f90", enzoFolder+fle+".F90")

		flist = ["opkda1", "opkda2", "opkdmain"]
		for fle in flist:
			shutil.move(buildFolder+fle+".f", enzoFolder+fle+".F")

		return
	############################################
	#prepare the patches if needed
	def patches(self):
		if(self.doFlash): self.flash_patch()
		if(self.doRamses): self.ramses_patch()
		if(self.doRamsesTH): self.ramsesTH_patch()
		if(self.doEnzo): self.enzo_patch()
		return

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
			#PATCHES DO NOT NEED MAKEFILE AND TEST.F90, and also if noExample is enabled
			if(not(self.doFlash or self.doRamses or self.doEnzo or self.doRamsesTH or self.noExample)):
				#TODO: add description in case of dust
				print "Call KROME in your code as:"
				if(useX):
					print "    call krome(x(:), gas_density, gas_temperature, time_step)"
				else:
					print "    call krome(x(:), gas_temperature, time_step)"
				print "where:" 
				print " x(:) is a real*8 array of size "+str(nmols)+(" of the mass fractions" if useX else\
					 " of number densities [1/cm3]")
				if(useX): print " gas_density  is the gas density in [g/cm3]"
				print " gas_temperature is the gas temperature in [K]"
				print " time_step is the integration time-step in [s]"
			
				if(not(self.isdry)):
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
