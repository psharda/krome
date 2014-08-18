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
# Others (alphabetically): F.A. Gianturco, T. Haugboelle, 
# J.Prieto, D.R.G. Schleicher, D. Seifried, E. Simoncini, 
# E. Tognelli
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
	useCoolingCompton = useCoolingExpansion = useShieldingDB96 = useShieldingWG11 = useCoolingCIE = useCoolingDISS = useCoolingFF = False
	#useCoolingZC = useCoolingZCp = useCoolingZSi = useCoolingZSip = useCoolingZO = useCoolingZOp = useCoolingZFe = useCoolingZFep = False
	useReverse = useCustomCoe = useODEConstant = cleanBuild = usePlainIsotopes = useDust = use_thermo = useStars = useNuclearMult = False
	usePhIoniz = useHeatingCompress = useHeatingPhoto = useHeatingChem = useDecoupled = useCoolingdH = useHeatingdH = useCoolingChem = False
	useHeatingCR = useHeatingPhotoAv = useHeatingPhotoDust = useHeatingXRay = useThermoToggle = False
	pedanticMakefile = useFakeOpacity = useConserve = useConserveE = noExample = useNLEQ = usePhotoOpacity = useXRay = False
	useX = has_plot = doIndent = useTlimits = useODEthermo = safe = doJacobian = True
	useDustGrowth = useDustSputter = useDustH2 = useDustT = checkThermochem = needLAPACK = useCoolCMBFloor = False
	doRamses = doRamsesTH = doFlash = doEnzo = wrapC = mergeTlimits = shortHead = isdry = useIERR = checkReverse = usePhotoInduced = False
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
	H2opacity = "" #H2 opacity model
	checkMode = "ALL" #conservation check mode (ALL | [CHARGE],[MASS]| NONE)
	RTOL = 1e-4 #default relative tolerance
	ATOL = 1e-20 #default absolute tolerance
	coolingQuench = -1e0 #if coolingQuench is negative cooling quench is not enabled, otherwise this is Tcrit
	dustArraySize = dustTypesSize = photoBins = 0
	maxord = 0
	dustTypes = []
	specs = []
	reacts = []
	constantList = []
	dummy = molec()
	coevars = dict() #variables in function coe() (krome_subs.f90)
	coevarsODE = dict() #variables in function fex() (krome_ode.f90)
	commonvars = [] #list of common variables
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
	nPhotoRea = 0 #number of photoreactions (for photobin array)
	full_cool = vars_cool = ""
	coolZ_functions = []
	coolZ_rates = []
	coolZ_vars_cool = []
	coolZ_poplevelvars = [] #population levels variables
	fcn_levs = [] #list of number of cooling levels found
	coolZ_nkrates = 0
	zcoolants = [] #list of cooling read from file (flag name, e.g CII)
	Zcools = [] #list of cooling read from file (species name, e.g. C+)
	anytabvars = [] #variable names for the tables
	anytabfiles = [] #file name for the tables
	anytabpaths = [] #paths for the tables
	anytabsizes = [] #sizes of the tables
	coolLevels = [] #levels employed for cooling, if empty uses all
	physVariables = [] #list of the phys variables (list of [variable_name, default_value_string])
	kModifier = [] #modifier lines that will be appended after the rate calculation
	odeModifier = [] #modifier lines that will be appended after the ODE calculation
	columnDensityMethod = "DEFAULT"
	ramses_offset = 3 #offset in the array for ramses
	coolFile = ["data/coolZ.dat"]
	fdbase = "data/database/"
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
		self.parser = argparse.ArgumentParser(description="KROME a package for astrochemistry and microphysics")
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
		self.parser.add_argument("-columnDensityMethod", metavar="method", help="use an alternative method to \
			N=1.8e21*(n*1e-3)**(2./3.) for column density calculation (N) from number density (n). Option available JEANS,\
			which employs Jeans length (l) as N=n*l.")
		self.parser.add_argument("-compressFluxes", action="store_true", help="in the ODE fluxes are stored in a single variable")
		self.parser.add_argument("-conserve", action="store_true", help="conserves the species total number and charge global\
			neutrality. Works with some limitations, please read the manual.")
		self.parser.add_argument("-conserveE", action="store_true", help="conserves the charge global neutrality only.")
		self.parser.add_argument("-coolFile", metavar='FILENAME', help="select the filename to be used to load external cooling. See\
			also tools/lamda2.py script for a LAMDA<->KROME converter. Default FILENAME is data/coolZ.dat, which contains\
			fine-strucutre atomic metal cooling for C,O,Si,Fe, and their first ions. It can also be a list of files comma-separated.")
		self.parser.add_argument("-cooling", metavar='TERMS', help="cooling options, TERMS can be ATOMIC, H2, HD, Z, DH, DUST, H2GP98,\
			COMPTON, EXPANSION, CIE, DISS, CI, CII, SiI, SiII, OI, OII, FeI, FeII, CHEM (e.g. -cooling=ATOMIC,CII,OI,FeI).\
			Note that further cooling options can be added when reading cooling function from file. If you want a complete list of\
			the available cooling options type -cooling=?")
		self.parser.add_argument("-coolLevels", metavar='MAXLEV', help="use only the levels up to MAXLEV (included), e.g. -coolLevels=3\
			Note that levels are zero-based (i.e. ground state is zero).")
		self.parser.add_argument("-coolingQuench", metavar='TCRIT', help="quenches the cooling when T<TCRIT with a tanh \
		 function.")
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
                        paritition function, ROT to keep into account the rotational partition function, EXACT to evaluate the\
                        adiabatic index accurately taking into account both contributions, or REDUCED to use only H2 and CO as diatomic\
			molecules (faster). Finally a custom F90 expression e.g. -gamma=\"1d0\"\
			can also be used. Default value is 5/3.",metavar="OPTION")
		self.parser.add_argument("-H2opacity", metavar="TYPE",help="use H2 opacity for H2 cooling, TYPE can be RIPAMONTI or OMUKAI")
		self.parser.add_argument("-heating", metavar='TERMS', help="heating options, TERMS can be COMPRESS, PHOTO, CHEM, DH, CR, PHOTOAV.\
			If you want a complete list of the available heating options type -heating=?")
		self.parser.add_argument("-ierr", action="store_true", help="same as -useIERR")
		self.parser.add_argument("-iRHS", action="store_true", help="implicit loop-based RHS (suggested for large systems)")
		self.parser.add_argument("-listAutomatics", action="store_true", help="list all the automatic reactions available")
		self.parser.add_argument("-maxord", help="max order of the BDF solver. Default (and maximum values) is 5.")
		self.parser.add_argument("-mergeTlimits", action="store_true", help="use the same reaction index for equivalent\
			reactions (same reactants and products) that have different temperature limits")
		self.parser.add_argument("-n", help="reaction network file", metavar='FILENAME')
		self.parser.add_argument("-network", help="same as -n", metavar='FILENAME')
		self.parser.add_argument("-nochargeCheck", action="store_true", help="skip reaction charge check")
		self.parser.add_argument("-noCheck", action="store_true", help="skip reaction charge and mass check. Equivalent to\
			-nomassCheck -nochargeCheck options.")
		self.parser.add_argument("-noExample", action="store_true", help="do not write test.f90 and Makefile in the build directory")
		self.parser.add_argument("-nomassCheck", action="store_true", help="skip reaction mass check")
		self.parser.add_argument("-noTlimits", action="store_true", help="ignore rate coefficient temperature limits.")
		self.parser.add_argument("-nuclearMult", action="store_true", help="keep into account reactants multeplicity, and modify\
			fluxes according to this. Intended for nuclear networks.")
		self.parser.add_argument("-options", metavar="filename", help="read the options from a file instead of command line\
			(in principle you can use both). See options_example file.")
		self.parser.add_argument("-pedantic", action="store_true", help="uses a pedantic Makefile (debug purposes)")
		self.parser.add_argument("-project", help="build everything in a folder called build_NAME instead of building all in the\
			default build folder. It also creates a NAME.kpj file with the krome input used.",metavar="NAME")
		self.parser.add_argument("-quote", action="store_true", help="print a citation and exit")
		self.parser.add_argument("-quotelist", action="store_true", help="print all the citations and exit")
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
		self.parser.add_argument("-photoBins", metavar="NBINS", help="define the number of frequency bins for the impinging radiation.")
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
		self.parser.add_argument("-useCoolCMBFloor", action="store_true", help="include a cooling floor given by the CMB temperature.\
			note that you must define Tcmb by using the subroutine krome_get_Tcmb(your_Tcmb) before calling krome.")
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
		self.parser.add_argument("-useIERR", action="store_true",help="use ierr in the interface with KROME to return errors instead\
			of stopping the exectution")
		self.parser.add_argument("-useN", action="store_true",help="use number densities (1/cm3) as input/ouput instead of fractions (#)")
		self.parser.add_argument("-useODEConstant", help="postpone an expression to each ODE. EXPRESSION must be a valid f90\
			expression (e.g. *3.d0 or +1.d-10)", metavar="EXPRESSION")
		self.parser.add_argument("-usePhIoniz", action="store_true", help="includes photochemistry (obsolete)")
		self.parser.add_argument("-usePhotoInduced", action="store_true", help="includes the photo-induced transitions in the calculation\
			of the cooling according to the choosen photon flux.")
		self.parser.add_argument("-usePhotoOpacity", action="store_true", help="computes photorates using opacity as a function of \
			the species densities and the photo cross sections, i.e. exp(-sum_i N_i*sigma_i). Column densities are computed\
			from density by using the local approximation N = 1.8e21*(n/1000)**(2/3) 1/cm2.")
		self.parser.add_argument("-usePlainIsotopes", action="store_true", help="use kA format for isotopes instead of [k]A format,\
			where k is the isotopic number and A is the atom name, e.g. krome looks for 14C instead of [14]C in the reactions file.")
		self.parser.add_argument("-useThermoToggle", action="store_true", help="include thermal calculation control. Use\
			krome_thermo_on and krome_thermo_off to switch on/off the thermal processes (i.e. cooling and heating). Default is on.")
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
			[argv.append(x) for x in ["-useN","-iRHS","-skipJacobian","-useCustomCoe=\"myCoe(:)\""]]
			filename = "networks/react_cloud"
		elif(args.test=="slowmanifold"):
			[argv.append(x) for x in ["-useN"]]
			filename = "networks/react_SM"
		elif(args.test=="auto"):
			[argv.append(x) for x in ["-photoBins=10","-useN"]]
			filename = "networks/react_auto"
		elif(args.test=="chianti"):
			[argv.append(x) for x in ["-photoBins=10","-useN","-useThermoToggle","-coolLevels=99999"]]
			[argv.append(x) for x in ["-cooling=CII,CIII,CIV,CV,CVI"]]
			[argv.append(x) for x in ["-coolFile=tools/coolChianti.dat"]]
			filename = "networks/react_chianti"
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
			[argv.append(x) for x in ["-cooling=ATOMIC,HD,H2", "-heating=PHOTO","-photoBins=10"]]
			filename = "networks/react_primordial_photoH2"
		elif(args.test=="collapse"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CONT,CHEM", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=FULL"]]
			filename = "networks/react_primordial3"
		elif(args.test=="collapseZ"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CI,CII,OI,OII,SiII,FeII,CONT,CHEM", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=FULL","-ATOL=1d-40","-maxord=1"]]
			filename = "networks/react_primordialZ2"
		elif(args.test=="collapseCO"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CI,CII,OI,CONT,CHEM", "-heating=COMPRESS,CHEM,CR,PHOTOAV,PHOTODUST"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=REDUCED","-ATOL=1d-10","-maxord=1","-useTabs"]]
			[argv.append(x) for x in ["-coolingQuench=10"]]
			filename = "networks/react_COthin"
		elif(args.test=="collapseZ_UV"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CI,CII,OI,OII,SiII,FeII,CONT,CHEM", "-heating=COMPRESS,CHEM,PHOTO"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=FULL","-photoBins=5","-usePhotoOpacity"]]
			filename = "networks/react_primordialZ2_UV"
		elif(args.test=="collapseZ_induced"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CI,CII,OI,OII,SiII,FeII,CONT,CHEM", "-heating=COMPRESS,CHEM,PHOTO"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=FULL","-photoBins=10","-usePhotoInduced"]]
			filename = "networks/react_primordialZ2"
		elif(args.test=="collapseUV"):
			[argv.append(x) for x in ["-cooling=H2,COMPTON,CIE,ATOMIC,FF", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-useN","-gamma=FULL"]]
			filename = "networks/react_primordial_UV"
                elif(args.test=="collapseUV_Xrays"):
			[argv.append(x) for x in ["-cooling=H2,CIE,DISS,HD,FF,COMPTON,CHEM", "-heating=COMPRESS,CHEM,XRAY"]]
			[argv.append(x) for x in ["-useN","-gamma=FULL","-shielding=WG11","-conserve","-H2opacity=OMUKAI"]]
			[argv.append(x) for x in ["-columnDensityMethod=JEANS"]]
			filename = "networks/react_xrays"
		elif(args.test=="collapseDUST"):
			[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE,DUST,HD", "-heating=COMPRESS,CHEM"]]
			[argv.append(x) for x in ["-H2opacity=RIPAMONTI","-useN","-gamma=FULL","-dust=1,C","-dustOptions=H2"]]
			filename = "networks/react_primordial"
                elif(args.test=="earlyUniverse"):
			[argv.append(x) for x in ["-cooling=H2GP98,COMPTON,EXPANSION"]]
			[argv.append(x) for x in ["-useN","-useFileIdx"]]
			filename = "networks/react_GP98"
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
		elif(args.test=="lamda"):
			[argv.append(x) for x in ["-useN","-coolFile=data/coolO2.dat", "-cooling=O2"]]
			[argv.append(x) for x in ["-useThermoToggle"]]
			filename = "networks/react_COthin"
		else:
			tests = ", ".join(sorted(os.walk('tests').next()[1]))
			print "ERROR: test \""+args.test+"\" not present!"
			print "Available tests are: "+tests
			sys.exit()

		self.filename = filename
		self.test_name = args.test

	##########################################
	def argparsing(self,argv):

		args = self.parser.parse_args() #return namespace from argv

		#list arguments if test
		if(args.test):
			print "This TEST is running with the following arguments:"
			for k in args.__dict__:
				arg = args.__dict__[k]
				if(arg): print " -"+k+" = "+str(arg)
			print " -n = "+self.filename
			print

		#use custom option file (load options from a file and append to argv)
		if(args.options):
			fopt = args.options.strip() #get filename
			print "Reading option -option="+fopt
			#check if option file exists
			if(not(file_exists(fopt))):
				print "ERROR: custom option file \""+fopt+"\" does not exist!"
				sys.exit()

			trues = ["T","TRUE","1","Y","YES","OK","YEP","SURE"]
			falses = ["F","FALSE","0","N","NO","KO","NOPE"]
			#read from file
			fho = open(fopt,"rb")
			for row in fho:
				srow = row.strip()
				if(srow==""): continue #skip blank lines
				if(srow[0]=="#"): continue #skip comments
				if(srow[:2]=="//"): continue #skip comments
				srow = srow.split("#")[0]
				srow = srow.split("//")[0]
				srow = srow.strip()
				#replace tabs
				srow = srow.replace("\t"," ")
				#replace double spaces
				while("  " in srow):
					srow = srow.replace("  "," ")
				if(srow[0]!="-"): srow = "-"+srow

				arow = srow.split()
				if(len(arow)==1): 
					sys.argv.append(arow[0].strip())
					continue
				elif(len(arow)==2):
					if(arow[1].strip().upper() in trues):
						sys.argv.append(arow[0].strip())
					elif(arow[1].strip().upper() in falses):
						continue
					else:
						sys.argv.append("=".join([x.strip() for x in arow]))
				else:
					print "ERROR: problems with option line in option file "+fopt
					print srow
					sys.exit()

			args = self.parser.parse_args() #return updated namespace

		#list all the automatic reactions available from the files in the fdbase folder and exit
		if(args.listAutomatics):
			os.path.isdir(self.fdbase)
			if(not(file_exists(self.fdbase))):
				print "ERROR: database directory "+self.fdbase+" not found!"
				sys.exit()
			file_list = [f for f in listdir(self.fdbase) if isfile(join(self.fdbase,f))]
			for fname in file_list:
				fname = self.fdbase+fname
				print "retriving reactions in "+fname
				fhauto = open(fname,"rb")
				icounta = 0
				reasa = prodsa = typea = ""
				for row in fhauto:
					srow = row.strip()
					if("@type:" in srow):
						reasa = prodsa = typea = ""
					if(reasa!="" and prodsa!="" and typea!=""):
						icounta += 1
						print str(icounta)+". ("+typea+") "+reasa+" -> "+prodsa
					if("@reacts:" in srow): 
						reasa = " + ".join([x.strip() for x in srow.replace("@reacts:","").split(",")])
					if("@prods:" in srow): 
						prodsa = " + ".join([x.strip() for x in srow.replace("@prods:","").split(",")])
					if("@type:" in srow): 
						typea = srow.replace("@type:","").strip()
				print
			sys.exit()
		
		#get a citation and exit
		if(args.quote):
			print "KROME is a quote random generator with some utility for astrochemistry."
			print "As requested a random citation:"
			get_quote()
			sys.exit()

		#get the list of the quotes and exit
		if(args.quotelist):
			print "KROME is a quote random generator with some utility for astrochemistry."
			print "As requested the complete list of the available citations:"
			get_quote(True)
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
		if(not(self.is_test) and args.network): self.filename = args.network
		#chech if reactions file exists
		if(args.n or self.is_test):
			if(not(os.path.isfile(self.filename))): die("ERROR: Reaction file \""+self.filename+"\" doesn't exist!")
		else:
			die("ERROR: you must define -n FILENAME or -network FILENAME, where FILENAME is the reaction file!")

		#read the coolFile
		if(args.coolFile):
			self.coolFile = args.coolFile.split(",")
			print "Reading option -coolFile (filename="+str(",".join(self.coolFile))+")"


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

		#method for column density calculation
		if(args.columnDensityMethod):
			allMethods = ["JEANS"]
			if(not(args.columnDensityMethod in allMethods)):
				sys.exit("ERROR: method for -columnDensityMethod must be one of "+(",".join(allMethods)))
			self.columnDensityMethod = args.columnDensityMethod

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
			#replace square brackets
			copydic = dict()
			for k,v in self.mass_dic.iteritems():
				copydic[k.replace("[","").replace("]","")] = v
			self.mass_dic = copydic
			self.atoms = [x.replace("[","").replace("]","") for x in self.atoms]
			
			
		#use photoionization from Verner et al. 1996 (no longer working)
		if(args.usePhIoniz):
			self.usePhIoniz = True
			print "Reading option -usePhIoniz (now obsolete, you can remove it)"

		#use photoionization 
		if(args.usePhotoOpacity):
			self.usePhotoOpacity = True
			print "Reading option -usePhotoOpacity (now obsolete, you can remove it)"

		#use cooling CMB floor 
		if(args.useCoolCMBFloor):
			self.useCoolCMBFloor = True
			if(not(args.cooling)):
				print "ERROR: option -useCoolCMBFloor needs at least one active cooling option. See -cooling="
				sys.exit()
			print "Reading option -useCoolCMBFloor"


		#use photo-induced cooling transitions 
		if(args.usePhotoInduced):
			self.usePhotoInduced = True
			if(not(args.photoBins)):
				print "ERRROR: -usePhotoInduced requires the option -photoBins=N enabled"
				print " where N is the number of photon bins employed."
				sys.exit()
			print "Reading option -usePhotoInduced"



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
		if(args.H2opacity):
			opacities = ["RIPAMONTI", "OMUKAI"]
			if(not(args.H2opacity in opacities)):
				print "ERROR: H2opacity must be one of "+(", ".join(opacities))+"."
				sys.exit()
			self.H2opacity = args.H2opacity.strip()
			print "Reading option -H2opacity="+self.H2opacity

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

		#include an if in the ODE for the thermal part
		if(args.useThermoToggle):
			self.useThermoToggle = True

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

		#do not write test.f90 and Makefile
		if(args.noExample):
			self.noExample = True
			print "Reading option -noExample"


		#set the number of photobins
		if(args.photoBins):
			self.photoBins = int(args.photoBins)
			self.usePhIoniz = True
			if(self.photoBins<0): die("ERRROR: number of frequency bins < 0!")
			print "Reading option -photoBins (NBINS="+str(self.photoBins)+")"

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
			allCools = ["ATOMIC","H2","HD","DH","DUST","FF","H2GP98","COMPTON","EXPANSION","CIE","CONT","CHEM","DISS","Z"]
			fileCools = [] #list of the cooling read from file
			#load additional coolings from file
			for fname in self.coolFile:
				partialFileCools = [] #for output purposes
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
						metal_name = srow.split(":")[1].strip()
						mol = parser(metal_name,self.mass_dic,self.atoms,self.thermodata)
						mname = mol.coolname
						if(mname in allCools):
							print "ERROR: conflict name for "+mname+", which is already present!"
							sys.exit()
						partialFileCools.append({"flag":mname,"name":metal_name})
						fileCools.append({"flag":mname,"name":metal_name}) #append to the list fo the available coolings
						allCools.append(mname) #append flag to the list of the coolants
				#write found coolants
				joinedCool = (", ".join([x["flag"] for x in partialFileCools]))
				if(len(partialFileCools)>0): 
					print "Cooling "+joinedCool+" available from "+fname
					allCools.append("FILE")

			if("?" in myCools):
				print "Available coolings are:", (", ".join(allCools))
				sys.exit()

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
			if("FF" in myCools): self.useCoolingFF = True
			if("DISS" in myCools): self.useCoolingDISS = True
			if("CONT" in myCools): self.useCoolingCont = True
			if("Z" in myCools): self.useCoolingZ = True

			#loop over metals loaded from file and search for them in the cooling flags provided by the user
			for met in fileCools:
				if((met["flag"] in myCools) or ("FILE" in myCools)):
					#print "Option "+met["flag"]+" will load data from "+fname
					self.useCoolingZ = True
					self.Zcools.append(met["name"]) #append metal name to the list of the requested species
					self.zcoolants.append(met["flag"]) #append cooling name to the list of the coolants
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

		if(args.coolLevels):
			self.coolLevels = [int(x) for x in range(int(args.coolLevels)+1)]
			if(int(args.coolLevels)<1): die("ERROR: coolLevels must be at least 1 (two levels)!")
			print "Reading option -coolLevels ("+str(self.coolLevels[0])+" to "+str(self.coolLevels[-1])+")"

		#cooling quenching
		if(args.coolingQuench):
			self.coolingQuench = format_double(args.coolingQuench)
			if(self.coolingQuench<0e0):
				die("ERROR: Tcrit for coolingQuench should be greater than zero!")
			print "Reading option -coolingQuench ("+str(self.coolingQuench)+")"
		
		#determine heating types
		if(args.heating):
			myHeat = args.heating.upper().split(",")
			myHeat = [x.strip() for x in myHeat]
			allHeats = ["COMPRESS","PHOTO","CHEM","DH","CR","PHOTOAV","PHOTODUST","XRAY"]
			for hea in myHeat:
				if(not(hea in allHeats)):
					die("ERROR: Heating \""+hea+"\" is unknown!\nAvailable heatings are: "+(", ".join(allHeats)))

			if("COMPRESS" in myHeat): self.useHeatingCompress = True
			if("PHOTO" in myHeat): self.useHeatingPhoto = True
			if("CHEM" in myHeat): self.useHeatingChem = True
			if("DH" in myHeat): self.useHeatingdH = True
			if("CR" in myHeat): self.useHeatingCR = True
			if("PHOTOAV" in myHeat): self.useHeatingPhotoAv = True
			if("PHOTODUST" in myHeat): self.useHeatingPhotoDust = True
			if("XRAY" in myHeat): self.useHeatingXRay = True

			self.use_thermo = True
			if(self.photoBins<=0 and self.useHeatingPhoto):
				print "ERROR: if you use photoheating you should include the number of photo-bins"
				print " by using the option -photoBins=NBINS"
				sys.exit()

			if("?" in myHeat):
				print "Available heatings are:", (", ".join(allHeats))
				sys.exit()


			print "Reading option -heating ("+(",".join(myHeat))+")"
	
		#force rwork size
		if(args.forceRWORK):
			myrwork = args.forceRWORK
			self.force_rwork = True
			print "Reading option -forceRWORK (RWORK="+str(myrwork)+")"

		#use custom function for coefficient instead of coe_tab()
		if(args.useCustomCoe):
			self.customCoeFunction = args.useCustomCoe.replace("\"","")
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
					#check if the current @var is allowed
					notAllowedVars = ["k","tgas","energy_ev"]
					for nav in notAllowedVars:
						if(nav.lower()==arow[0].lower()):
							sys.exit("ERROR: you can't use "+nav+" as an @var variable")
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
						print "ERROR: NASA polynomials!"
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
			'Ne':10.*(menp),
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

		#excited levels of some molecules. add here if needed
		for i in range(9):
			mass_dic['CH2_'+str(i+1)] = 6.*(menp) + 2.*(me+mp)
			mass_dic['SO2_'+str(i+1)] = 16.*(menp) + 2.*8.*(menp)

		#build isotopes (including some non-esistent) as [n]A
		# with -usePlainIsotopes build as nA
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

		#prepare mass dictionary
		self.mass_dic = dict([[k.upper(),v] for (k,v) in mass_dic.iteritems()])
		#sort dictionary, longest first. note that even if it is called 
		# atoms, this contains also other chemical formula parts, as GRAIN, PAH, and so on... 
		self.atoms = sorted(mass_dic, key = lambda x: len(x),reverse=True)

	
	#################################################
	#read the reaction file 
	def read_file(self):
		skipDup = self.skipDup
		filename = self.filename
		atoms = self.atoms
		mass_dic = self.mass_dic
		thermodata = self.thermodata
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
		specs = [] #list of species as mol objects
		reacts = [] #list of reactions as react objects
		reags = [] #list of reagents for already found
		prods = [] #list of prods for already found
		idxs = [] #list of index for already found
		noTabNextBlock = False #default for blocks of reactions
		inCRblock = False #block of CR reactions
		inPhotoBlock = False #block of photo reactions with xsection function
		inXRayBlock = False #block of xray reactions
		inReactionModifierBlock = False #block of modifier expression for the computed coefficients
		inOdeModifierBlock = False #block of modifier expression for the ODE dn/dT
		noTabBlockStored = noTabNextBlock #store the noTabNextBlock array before inPhotoBlock to restore it

		#read the size of the file in lines (skip blank and comments)
		# to have a rough idea of the size
		fh = open(filename,"rb")
		line_count = 0
		allrows = []
		for row in fh:
			if(row.strip()==""): continue
			if(row.strip()[0]=="#"): continue
			line_count += 1
			allrows.append(row.strip())
		fh.close()
	
		#warning if the number of lines exceed a certain limit
		if(line_count>1000): print "Found "+str(line_count)+" lines! It takes a while..."

		fsh_found = False #search for fsh variable for shielding if needed
		#start reading file stored in the loop above
		isComment = False #flag for comment block
		noTabNext = False #flag for use tabs for the next reaction 
		for row in allrows:
			srow = row.strip() #stripped row
			if(srow.strip()==""): continue #looks for blank line
			if(srow[0]=="#"): continue #looks for comment line
			if(len(srow)>1):
				if(srow[0:2]=="//"): continue #looks for comment line
				if(srow[0:2]=="/*"): isComment = True #start multiline comment

			#end multiline comment
			if("*/" in srow):
				isComment = False
				continue
			if(isComment): continue #skip if in comment block


			#search for final expression to modify the coefficients (stop)
			if(srow.lower()=="@reactionmodifier_stop" or srow.lower()=="@reactionmodifier_end"):
				inReactionModifierBlock = False
				continue #SKIP (not a reaction)

			#search for final expression to modify the ODE (stop)
			if(srow.lower()=="@odemodifier_stop" or srow.lower()=="@odemodifier_end"):
				inOdeModifierBlock = False
				continue #SKIP (not a reaction)

			#store coefficient modifier
			if(inReactionModifierBlock):
				if("@" in srow):
					print srow
					sys.exit("ERROR: @ expressions are not allowed inside the @reactionModifier block!")
				self.kModifier.append(srow)
				continue #skip: modifier line is not a reaction

			#store ode modifier
			if(inOdeModifierBlock):
				if("@" in srow):
					print srow
					sys.exit("ERROR: @ expressions are not allowed inside the @ODEModifier block!")
				self.odeModifier.append(srow)
				continue #skip: modifier line is not a reaction


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
				#check if the current @var is allowed
				notAllowedVars = ["k","tgas","energy_ev"]
				for nav in notAllowedVars:
					if(nav.lower()==arow[0].lower()):
						sys.exit("ERROR: you can't use "+nav+" as an @var variable")

				if(arow[0] in self.coevars): continue #skip already found variables
				self.coevars[arow[0]] = [ivarcoe,arow[1]]
				ivarcoe += 1 #count variables to sort
				continue #SKIP: a variable line is not a reaction line

			#search for common variables
			if("@common:" in srow):
				arow = srow.replace("@common:","").split(",")
				for x in arow:
					commonvar = x.strip()
					if(commonvar.split("_")[0].lower()!="user"):
						print "ERROR: to avoid conflicts common variables with @common should begin with user_"
						print " you provided: "+commonvar
						print " it should be: user_"+commonvar
						sys.exit()
					if(commonvar in self.commonvars): continue #skip if already present
					self.commonvars.append(commonvar) #add to the global array
				continue #skip: a common is not a reaction line

			#search for ghost species
			if("@ghost:" in srow):
				ghost = srow.replace("@ghost:","").strip()
				aghost = ghost.split(",")
				for ghost in aghost:
					print "Found ghost species: "+ghost
					mol = parser(ghost,mass_dic,atoms,self.thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
				continue #SKIP: a ghost line is not a reaction line

			#search for table pragma
			if("@tabvar:" in srow):
				atab = srow.replace("@tabvar:","").split("=")
				if(len(atab)!=2):
					print "ERROR: wrong format, it should be @tabvar:varname=file,var1,var2"
					print " You provided: "+srow.strip()
					sys.exit()
				aatab = atab[1].split(",")
				if(len(aatab)!=3):
					print "ERROR: wrong format, it should be @tabvar:varname=file,var1,var2"
					print " You provided: "+srow.strip()
					sys.exit()
				
				mytabvar = atab[0].strip()
				mytabpath = aatab[0].strip().replace("\"","")
				mytabxxyy = aatab[1]+","+aatab[2]
				#updates anytab arrays
				create_tabvar(mytabvar,mytabpath,mytabxxyy,self.anytabvars,self.anytabfiles,self.anytabpaths,\
					self.anytabsizes,self.coevars,ivarcoe)
				continue #this is not a reaction line
				
			#search for format string
			if("@format:" in srow):
				idxFound = tminFound = tmaxFound = rateFound = qeffFound = False
				hasFormat = True #format flag
				srow = srow.replace("@format:","") #remove 
				#print "Found custom format: "+srow
				arow = srow.split(",") #split format line
				#check format (at least 6 elements)
				if(len(arow)<5):
					print "ERROR: format line must contains at least 5 elements"
					print " idx,R,P,P,rate"
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
			if(srow.lower()=="@notabnext" or srow.lower()=="@notab_next"):
				noTabNext = True
				continue #SKIP (not a reaction)
			#if requested the next reaction will not uses tabs for a BLOCK of reactions
			if(srow.lower()=="@notab_begin" or srow.lower()=="@notab_start"):
				noTabNext = noTabNextBlock = True
				continue #SKIP (not a reaction)
			#if requested the next reaction will not uses tabs for a BLOCK of reactions
			if(srow.lower()=="@notab_end" or srow.lower()=="@notab_stop"):
				noTabNext = noTabNextBlock = False
				continue #SKIP (not a reaction)

			#start a CR reaction block
			if(srow.lower()=="@cr_start" or srow.lower()=="@cr_begin"):
				inCRblock = True
				noTabBlockStored = noTabNextBlock
				noTabNext = noTabNextBlock = True
				continue #SKIP (not a reaction)

			#start a CR reaction block
			if(srow.lower()=="@cr_stop" or srow.lower()=="@cr_end"):
				inCRblock = False
				noTabNext = noTabNextBlock = noTabBlockStored #restore the noTabNextBlock value before entering CR block
				continue #SKIP (not a reaction)

			#start a photo reaction block
			if(srow.lower()=="@photo_start" or srow.lower()=="@photo_begin"):
				if(self.photoBins<=0):
					print "ERROR: you are using "+srow.lower()+" in your reaction file"
					print " with zero photo-bins. Use -photoBins=NBINS option."
					sys.exit()
				inPhotoBlock = True
				noTabBlockStored = noTabNextBlock
				noTabNext = noTabNextBlock = True
				continue #SKIP (not a reaction)

			#start a photo reaction block
			if(srow.lower()=="@photo_stop" or srow.lower()=="@photo_end"):
				inPhotoBlock = False
				noTabNext = noTabNextBlock = noTabBlockStored #restore the noTabNextBlock value before entering inPhotoBlock
				continue #SKIP (not a reaction)

			#start an XRAY reaction block
			if(srow.lower()=="@xray_start" or srow.lower()=="@xray_begin"):
				inXRayBlock = True
				self.useXRay = True
				noTabBlockStored = noTabNextBlock
				noTabNext = noTabNextBlock = True
				continue #SKIP (not a reaction)

			#start an XRAY reaction block
			if(srow.lower()=="@xray_stop" or srow.lower()=="@xray_end"):
				inXRayBlock = False
				noTabNext = noTabNextBlock = noTabBlockStored #restore the noTabNextBlock value before entering inXRayBlock
				continue #SKIP (not a reaction)

			#search for final expression to modify the coefficients (start)
			if(srow.lower()=="@reactionmodifier_start" or srow.lower()=="@reactionmodifier_begin"):
				inReactionModifierBlock = True
				continue #SKIP (not a reaction)

			#search for final expression to modify the coefficients (start)
			if(srow.lower()=="@odemodifier_start" or srow.lower()=="@odemodifier_begin"):
				inOdeModifierBlock = True
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
			rcount += 1 #count the total number of reaction found
			printStep = int(line_count/10)
			if((line_count>1000) and (rcount%printStep==0)): print str(round(rcount*1e2/line_count))+"%"

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
			if("krome_fshield" in myrea.krate.lower()): fsh_found = True
			
			if(qeffFound): myrea.qeff = arow[iqeff]


			#if(self.useCustomCoe): myrea.krate = "0.d0" #when custom function is used standard coefficient are set to zero
			#loop over reactants to grep molecules
			for r in reactants:
				if(r.strip()=="G" and not(self.use_photons)): continue
				if(r.strip()=="E-"): r = "E"
				if(r.strip()!=""):
					mol = parser(r,mass_dic,atoms,thermodata)
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
					mol = parser(p,mass_dic,atoms,thermodata)
					if(not(mol.name in spec_names)):
						spec_names.append(mol.name)
						specs.append(mol)
					mol.idx = spec_names.index(mol.name) + 1
					myrea.products.append(mol) #add molecule object to products

			#increases the index of the photoreaction
			if(inPhotoBlock):
				self.nPhotoRea += 1
				myrea.idxph = self.nPhotoRea

			myrea.build_verbatim() #build reaction as string (e.g. A+B->C)
			#myrea.reactants = sorted(myrea.reactants, key=lambda r:r.idx) #sort reactants
			#myrea.products = sorted(myrea.products, key=lambda p:p.idx) #sort products
			myrea.build_RHS(self.useNuclearMult) #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
			myrea.build_phrate(inPhotoBlock) #build photoionization rate
			myrea.check(self.checkMode) #check mass and charge conservation
			myrea.group = group #add the group to the reaction
			myrea.canUseTabs = not(noTabNext) #check if this reaction can use tabs or not
			if(myrea.krate.count("(")!=myrea.krate.count(")")):
				print "ERROR: unbalanced brackets in reaction "+str(myrea.idx)
				print " "+myrea.verbatim
				print " rate = "+myrea.krate
				print " this is the corresponding line in the reaction file"
				print srow
				sys.exit()

			myrea.isCR = inCRblock #is a CR reaction
			myrea.isXRay = inXRayBlock #is an XRAY reaction

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
			del myrea,row
			if(not(noTabNextBlock)): noTabNext = False #return to default value when outside a block
	
		if((self.useShieldingDB96 or self.useShieldingWG11) and not(fsh_found)):
			print
			print "WARNING: no krome_fshield(n(:),Tgas) variable found in rate coefficient"
			print " even if shielding option is enabled."
			print " Please check your network file!"
			a = raw_input("Any key to continue q to quit... ")
			if(a=="q"): print sys.exit()

		if(noTabNextBlock): 
			print "ERROR: block of skipped reaction still open!"
			print "Add @noTab_stop or @noTab_end"
			sys.exit()

		if(skipDup): 
			fdup.close()
			print "Skipped duplicated reactions:",skipped_dupl

		#check file format
		if(not(found_one)):
			die("ERROR: no valid reactions found in file \""+filename+"\"")
		if(unmatch_idx):
			print "WARNING: index in \""+filename+"\" are not sequential!"

		#prepares xray rates including self-shielding and secondary process
		xrayHFound = xrayHeFound = False
		for x in reacts:
			if(not(x.isXRay)): continue
			if(not("auto" in x.krate)): continue
			ivarcoe = len(self.coevars)
			fake_ivarcoe = 0
			fake_coevars = dict()
			addVarCoe("ncolH","num2col(n(idx_H),n(:))",self.coevars)
			addVarCoe("ncolHe","num2col(n(idx_He),n(:))",self.coevars)
			addVarCoe("logHe","log10(ncolHe)",self.coevars)
			addVarCoe("logH","log10(ncolH)",self.coevars)
			addVarCoe("xe","n(idx_e) / (get_Hnuclei(n(:)) + 1d-40)",self.coevars)

			#updates anytab arrays
			if(x.reactants[0].name=="H"):
				mytabvar = "user_xray_H"
				mytabpath = "data/ratexH.dat"
				mytabxxyy = "logH,logHe-logH"

				create_tabvar(mytabvar,mytabpath,mytabxxyy,self.anytabvars,self.anytabfiles,self.anytabpaths,\
					self.anytabsizes,self.coevars)

				addVarCoe("phiH",".3908d0*(1e0-xe**.4092)**1.7592 * 327.832286034056d0",self.coevars)
				addVarCoe("ratexH"," 1d1**user_xray_H",self.coevars)

				xrayHFound = True
				autoRateXray = "ratexH * (1d0+phiH) + n(idx_He)/(n(idx_H)+1d-40) * ratexHe * phiH"
				x.krate = autoRateXray + "* J21xray"
				print "H xray ionization found!"

				#heating tabs H
				mytabvar = "user_xheat_H"
				mytabpath = "data/heatxH.dat"
				mytabxxyy = "logH,logHe-logH"
				create_tabvar(mytabvar,mytabpath,mytabxxyy,self.anytabvars,self.anytabfiles,self.anytabpaths,\
					self.anytabsizes,fake_coevars)

			elif(x.reactants[0].name.lower()=="he"):
				mytabvar = "user_xray_He"
				mytabpath = "data/ratexHe.dat"
				mytabxxyy = "logH,logHe-logH"

				create_tabvar(mytabvar,mytabpath,mytabxxyy,self.anytabvars,self.anytabfiles,self.anytabpaths,\
					self.anytabsizes,self.coevars)

				addVarCoe("phiHe",".0554d0*(1e0-xe**.4614)**1.666 * 180.793458763612d0",self.coevars)
				addVarCoe("ratexHe"," 1d1**user_xray_He",self.coevars)

				autoRateXRay = "ratexHe * (1d0+phiHe) + n(idx_H)/(n(idx_He)+1d-40) * ratexH * phiHe"
				x.krate = autoRateXray + "* J21xray"
				xrayHeFound = True
				print "He xray ionization found!"

				#heating tabs He
				mytabvar = "user_xheat_He"
				mytabpath = "data/heatxHe.dat"
				mytabxxyy = "logH,logHe-logH"
				create_tabvar(mytabvar,mytabpath,mytabxxyy,self.anytabvars,self.anytabfiles,self.anytabpaths,\
					self.anytabsizes,fake_coevars)

			else:
				print "ERROR: xray reaction not tabulated!"
				print " "+x.verbatim
				print " remove it from the chemical network or provide non-automatic rate."
				print " Note that you should also provide the heating tab if needed."
				sys.exit()

		#check if both (H and He) xray reactions are found, since tables are H and He dependant
		if(xrayHeFound!=xrayHFound):
			print "ERROR: for xrays you must include both H and He reaction in your reaction files, e.g.:"
			print " @format:idx,R,P,P,Tmin,Tmax,rate"
			print " 7,H,H+,E,NONE,NONE,auto"
			print " 8,He,He+,E,NONE,NONE,auto"
			print " where indexes are arbitrary"
			sys.exit()

		#check if both (H and He) xray reactions are found, when xray heating is enabled
		if((not(xrayHeFound) and not(xrayHFound)) and self.useHeatingXRay):
			print "ERROR: with XRAY heating option you must include both H and He reaction"
			print " in your reaction files, e.g.:"
			print " @format:idx,R,P,P,Tmin,Tmax,rate"
			print " 7,H,H+,E,NONE,NONE,auto"
			print " 8,He,He+,E,NONE,NONE,auto"
			print " where indexes are arbitrary"
			sys.exit()

		#check for automatic reactions
		autoFound = False
		for rea in reacts:
			if(rea.kphrate!=None):
				if(rea.kphrate.lower().strip()=="auto"):
					autoFound = True
					break
			if(rea.krate!=None):
				if(rea.krate.lower().strip()=="auto"):
					autoFound = True
					break
		
		#search auto reaction in the database
		if(autoFound):
			autoreacts = [] #dbase array contains dictionary with reaction data
			fdbase = self.fdbase
			print "Automatic reactions found, searching in "+fdbase
			if(not(os.path.isdir(fdbase))):
				print "ERROR: folder "+fdbase+" not found!"
				sys.exit()
			file_list = [f for f in listdir(self.fdbase) if isfile(join(self.fdbase,f))]
			for fname in file_list:
				fname = fdbase + fname
				fhdbase = open(fname,"rb") #open the database
				#load the database into an array of dictionaries
				for row in fhdbase:
					srow = row.strip()
					if(srow==""): continue #skip blank
					if(srow[0]=="#"): continue #skip comments
					#each reaction block starts with @type, init the reaction dictionary
					if("@type:" in srow):
						myrea = dict()
					myrea.update(at_extract(srow)) #append to the dictionary
					#each reaction block ends with @rate, append to the main database array
					if("@rate:" in srow):
						autoreacts.append(myrea)
			#loop on the reactions to find auto
			for i in range(len(reacts)):
				rea = reacts[i]
				if(rea.kphrate==None):
					if(rea.krate.lower().strip()!="auto"): continue
				else:
					if(rea.krate.lower().strip()!="auto" and rea.kphrate.lower().strip()!="auto"): continue
				dbFound = False
				#loop on autoreactions
				for autorea in autoreacts:
					autop = [x.upper().strip() for x in autorea["prods"].split(",")] #list of prods
					autor = [x.upper().strip() for x in autorea["reacts"].split(",")] #list of reacts
					if(sorted([x.name for x in rea.reactants])!=sorted(autor)): continue
					if(sorted([x.name for x in rea.products])!=sorted(autop)): continue
					dbFound = True
					#handle photochemistry
					if(rea.kphrate=="auto"):
						reacts[i].kphrate = autorea["rate"]
					else:
						reacts[i].krate = autorea["rate"]

					reacts[i].Tmin = autorea["limits"].split(",")[0].strip()
					reacts[i].Tmax = autorea["limits"].split(",")[1].strip()
					print "automatic reaction found!",rea.verbatim
					break
				#error if automatic reaction not found
				if(not(dbFound)):
					print "ERROR: reaction not found in the automatc database!"
					print rea.verbatim,[x.name for x in rea.reactants]
					print "you can:"
					print "1. remove it from your network"
					print "2. provide a non-automatic reaction rate"
					print "3. add to the databse "+fdbase
					sys.exit()


		#count reactions with unique index
		idxs = []
		nrea = 0
		for rea in reacts:
			if(rea.idx in idxs):
				print "skipping "+rea.verbatim+" since index already used!"
				continue #skip reactions same index
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
	#define the phys_ variables (will be used in krome_commons and 
	# in krome_user to create the get and set functions) 
	def definePhysVariables(self):
		#variables are list [name, default_value_string]
		#note that phys_ will be prepended 
		self.physVariables = [["Tcmb", "2.73d0"],
			["zredshift", "0d0"],
			["orthoParaRatio", "3d0"]]
	
	#####################################################
	def photo_warnings(self):
		if(self.is_test): return #skip warning if test mode
		if(self.usePhIoniz or self.useHeatingPhoto):
			print "************************************************"
			print "REMINDER: note that, since you are using photon-based"
			print " options, you need to initialize the machinery from"
			print " your main file! Read the manual for further details."
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
		if(self.nPhotoRea>0): print "Photo reactions found: ",self.nPhotoRea


	########################
	def uniq(self,a):
		u = []
		for x in a:
			if(not(x in u)): u.append(x)
		return u

	###############################################
	def dumpNetwork(self):

		buildFolder = self.buildFolder
		#create build folder if not exists
		#(this block is also in prepareBuild method)
		if(not(os.path.exists(buildFolder))):
			os.mkdir(buildFolder)
			print "Created "+buildFolder
		
		#dump species to log file
		fout = open(self.buildFolder+"species.log","w")
		fout.write("#This file contains a list of the species used with their indexes\n")
		fout.write("\n")
		idx = 0
		for mol in self.specs:
			idx += 1
			fout.write(str(idx)+"\t"+mol.name+"\t"+mol.fidx+"\n")
		fout.close()
		print "Species list saved in "+self.buildFolder+"species.log"

		#dump species to gnuplot initialization
		fout = open(self.buildFolder+"species.gps","w")
		fout.write("#This file is a script to initialize the species index in krome gnuplot\n")
		idx = 0
		fout.write("\n")
		fout.write("\n")
		fout.write("if(!exists(\"nkrome\")) print \"ERROR: first you must set the value of the offset nkrome\"\n")
		inits = []
		for mol in self.specs:
			idx += 1
			inits.append("krome_"+mol.fidx+" = "+str(idx)+" + nkrome")
		fout.write(("\n".join(inits))+"\n")
		fout.write("print \"All variables set as e.g. krome_idx_H2\"\n")
		fout.write("print \"plot 'your_file' u 1:(column(krome_idx_H2))\"\n")
		fout.write("print \" the offset is nkrome=\",nkrome\n")
		fout.close()
		print "Species index initialization for gnuplot in "+self.buildFolder+"species.gps"

		#dump heating and cooling index initialization for gnuplot
		fout = open(self.buildFolder+"heatcool.gps","w")
		fout.write("#This file is a script to initialize the heating and cooling index in krome gnuplot\n")
		idx = 0
		fout.write("\n")
		fout.write("\n")
		fout.write("if(!exists(\"nkrome_heatcool\")) print \"ERROR: first you must set the value of the offset nkrome_heatcool\"\n")
		idxcools = get_cooling_index_list()
		idxheats = get_heating_index_list()
		inits = []
		for idx in idxcools:
			inits.append("krome_"+idx+" + nkrome_heatcool")
		for idx in idxheats:
			inits.append("krome_"+idx+" + nkrome_heatcool")
		fout.write(("\n".join(inits))+"\n")
		fout.write("print \"All variables set as e.g. krome_idx_cool_H2\"\n")
		fout.write("print \"plot 'your_file' u 1:(column(krome_idx_cool_H2))\"\n")
		fout.write("print \" the offset is nkrome_heatcool=\",nkrome_heatcool\n")
		fout.close()
		print "Heating cooling index init for gnuplot in "+self.buildFolder+"heatcool.gps"

	
		#dump reactions to log file
		
		fout = open(self.buildFolder+"reactions.log","w")
		idx = maxprod = maxreag = 0
		for rea in self.reacts:
			idx += 1
			rcount = pcount = 0
			#search for the maximum number of reactants and products
			for r in rea.reactants:
				if(r.name!=""): rcount += 1
			maxreag = max(maxreag,rcount)
			for p in rea.products:
				if(p.name!=""): pcount += 1
			maxprod = max(maxprod,pcount)
			fout.write(str(rea.idx)+"\t"+rea.verbatim+"\n")
		fout.close()

		
		#dump network to dot file
		fout = open(self.buildFolder+"network.dot","w")
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
		print "Reactions saved in "+self.buildFolder+"reactions.log"
	
	##############################################
	#write the C header if needed
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
	#show info about ODE system
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

	###################################
	def convertMetal2Roman(self,arg_metal):
		if("+" in arg_metal):
			return arg_metal.replace("+","") + int_to_roman(arg_metal.count("+")+1)
		elif("-" in arg_metal):
			return arg_metal.replace("-","") + "m"+int_to_roman(arg_metal.count("-")+1)
		else:
			return arg_metal+"I" #is not an ion

	###################################
	def convertMetal2F90(self,arg_metal):
		if("+" in arg_metal):
			return arg_metal.replace("+","j")
		elif("-" in arg_metal):
			return arg_metal.replace("-","k")
		else:
			return arg_metal #is not an ion
	
		
	#####################################
	#alternative for reading cooling data from file
	def createZcooling(self):
		cooling_data = dict() #structure with all the data
		index_count = 0 #1-based index for all the rates of the cooling
		TgasUpperLimit = TgasLowerLimit = "" #temperature limits for ALL the reactions (empty=no limits)
		#PART1: read the data from file
		for fname in self.coolFile:
			#read the file fname
			skip_metal = False #flag to skip metals not included in the self.zcoolants list
			print "******************"
			print "Reading coolants from "+fname+"..."
			fh = open(fname,"rb") #open file
			#loop on file lines
			for row in fh:
				srow = row.strip()
				#skip blank lines and comments
				if(srow==""): continue
				if(srow[0]=="#"): continue
				if(srow[:2]=="//"): continue
				#read metal name from file and init data structures
				if("metal:" in srow):
					srow = srow.replace("metal:","").strip()
					mol = parser(srow,self.mass_dic, self.atoms, self.thermodata)
					metal_name = mol.coolname
					metal_name_f90 = self.convertMetal2F90(srow) #metal name in f90 format (e.g. C+ as Cj)
					levels_data = dict() #data of the levels (energy, degeneration), key=level number
					trans_data = dict() #transition data (up, down, Aij), key=up->down (e.g. "4->2")
					rate_data = dict() #rate data (up,dowm,rate,collider), key=up->down_coolider (e.g. "4->2_H")
					skip_metal = False
					needOrthoPara = False #need extra variables for ortho/para H2
					#check if the metal name is in the self.zcoolant list to skip it
					if(not(metal_name.upper() in [x.upper() for x in self.zcoolants])): skip_metal = True
					continue
				if(skip_metal): continue #skip the metal if necessary (not included in the self.zcoolant list)

				#if level pragma found read level data
				if("level" in srow):
					srow = srow.replace("level","").replace(":",",") #use only comma to separate and remove "level"
					arow = [x.replace("d","e") for x in srow.split(",")] #split using comma and floating format conversion

					#skip if levels are not requested (-coolLevels)
					cond1 = (len(self.coolLevels)!=0) #condition1: the list of the requested level should be not empty
					cond2 = not(int(arow[0]) in self.coolLevels) #condition2: the level should be in the list
					if(cond1 and cond2): continue #skip if levels are not listed in requested levels

					#add level data to the dictionary
					levels_data[int(arow[0])] = {"energy":float(arow[1]), "g":float(arow[2])}
					continue
				
				#if 3 elemets is transistion data
				if("->" in srow):
					srow = srow.replace("->",",") #replace the arrow with a comma (easier to split)
					arow = [x.replace("d","e") for x in srow.split(",")] #split using comma and floating format conversion
					trans_name = arow[0].strip()+"->"+arow[1].strip() #prepare the key as up->down (e.g. "4->1")
					#skip if levels are not requested (-coolLevels)
					cond1 = (len(self.coolLevels)!=0) #condition1: the list of the requested level should be not empty
					cond2 = not(int(arow[0]) in self.coolLevels) #condition2: the up level should be in the list
					cond3 = not(int(arow[1]) in self.coolLevels) #condition3: the low level should be in the list
					if(cond1 and (cond2 or cond3)): continue #skip if levels are not listed in requested levels

					#store data
					kboltzmann_eV = 8.617332478e-5 #eV/K
					hplanck_eV = 4.135667516e-15 #eV*s
					clight =  2.99792458e10 #cm/s
					de_eV = kboltzmann_eV * abs(levels_data[int(arow[0])]["energy"] - levels_data[int(arow[1])]["energy"])
					#if wavelength in angstrom is available and positive use it to compute delta energy
					if(len(arow)>3):
						wvl = float(arow[3])/1e8 #cm
						de_eV = 0e0
						if(wvl>0e0):
							de_eV = hplanck_eV*clight/wvl #eV

					#compute pre-factor for photo-induced transitions
					if(de_eV!=0e0):
						preB = .5e0*(hplanck_eV*clight)**2/(de_eV)**3 #cm2/eV
					else:
						preB = 0e0

					#compute gi and gj
					gi = levels_data[(int(arow[0]))]["g"]
					gj = levels_data[(int(arow[1]))]["g"]

					#append transition dictionary to transition data
					trans_data[trans_name] = {"up":int(arow[0]), "down":int(arow[1]), "Aij":float(arow[2]),\
						"Bij":float(arow[2])*preB, "denergy_eV":de_eV, "denergy_K":de_eV/kboltzmann_eV,\
						"Bji":float(arow[2])*preB/gj*gi}
					continue
		
				#if more than 3 elements is a rate data
				if(len(srow.split(","))>3):
					arow = srow.split(",")
					arow[3] = (",".join(arow[3:])) #join the last element to cope with comma separated f90 functions
					#key for the transition as up->down_collider (e.g. "6->1_H")
					trans_name = arow[1].strip()+"->"+arow[2].strip()+"_"+arow[0].strip()

					#skip if levels are not requested (-coolLevels)
					cond1 = (len(self.coolLevels)!=0) #condition1: the list of the requested level should be not empty
					cond2 = not(int(arow[1]) in self.coolLevels) #condition2: the up level should be in the list
					cond3 = not(int(arow[2]) in self.coolLevels) #condition3: the low level should be in the list
					if(cond1 and (cond2 or cond3)): continue #skip if levels are not listed in requested levels

					#increase the number of the reactions found
					index_count += 1

					#store the reaction rate in the list
					rate_comment = "!"+arow[1].strip()+"->"+arow[2].strip()+", "+metal_name+" - "+arow[0].strip()+"\n"
					self.coolZ_rates.append(rate_comment+"k("+str(index_count)+") = "+arow[3].strip())

					#store data
					rate_data[trans_name] = {"up":int(arow[1]), "down":int(arow[2]), "collider":arow[0].strip(),\
						"rate":index_count}

					#check if para/ortho needed
					if(arow[0].strip()=="H2or" or arow[0].strip()=="H2pa"): needOrthoPara = True

					continue

				#check for if conditions on the rate (note that index_count is not incremented)
				if(srow[:3].lower()=="if("):
					arow = srow.split(":")
					self.coolZ_rates.append(arow[0].strip()+" k("+str(index_count)+") = "+arow[1].strip())
					continue

				#read @var
				if("@var:" in srow):
					srow = srow.replace("@var:","").strip()
					#read extra variables as [variable_name, expression]
					self.coolZ_vars_cool.append([x.strip() for x in srow.split("=")])
					continue

				#look for @TgasUpperLimit
				if("@tgasupperlimit:" in srow.lower()):
					TgasUpperLimit = srow.split(":")[1].strip()
					#append the upper limit for the temperature as the first reaction rate
					self.coolZ_rates.append("if(Tgas.ge."+TgasUpperLimit+") return\n")

				#end of data, hence store
				if(("endmetal" in srow) or ("end metal" in srow)):
					#store all the data for the given metal_name
					cooling_data[metal_name] = {"levels_data":levels_data, "trans_data":trans_data, "rate_data":rate_data,\
						"metal_name_f90":metal_name_f90, "needOrthoPara":needOrthoPara}
					#pprint(cooling_data)


		#PART2: use data to prepare cooling routine
		#prepare the functions for the cooling looping on metals (which are the key of the cooling_data dictionary)
		for cur_metal,cool_data in cooling_data.iteritems():
			metal_name = cur_metal #alias for metal name
			metal_name_f90 = cool_data["metal_name_f90"] #name in f90 style
			level_list = cool_data["levels_data"].keys() #store the list of the levels as integer values (e.g. [0,1,3])
			#make a local copy of the dictionaries (easy to handle)
			levels_data = cooling_data[cur_metal]["levels_data"]
			trans_data = cooling_data[cur_metal]["trans_data"]
			rate_data = cooling_data[cur_metal]["rate_data"]
			rate_data_rev = dict() #init a dictionary to store the inverse reactions
			print "Prepearing "+cur_metal+" (levels: "+str(len(levels_data))+", transitions: "+str(len(rate_data))+")"

			#PART2.1: prepare excitation rates from de-excitations
			#loop on rates to prepare the reverse (excitation rates)
			for trans_name, r_data in rate_data.iteritems():
				#inverse transition name
				trans_name_rev = str(r_data["down"])+"->"+str(r_data["up"])+"_"+r_data["collider"]
				g_up = levels_data[r_data["up"]]["g"] #upper level
				g_down = levels_data[r_data["down"]]["g"] #lower level
				#delta energy Ej-Ei
				deltaE = abs(levels_data[r_data["up"]]["energy"] - levels_data[r_data["down"]]["energy"])
				deltaE = ("%e" % deltaE).replace("e","d") #f90ish format for deltaE
				#increase the number of the reactions found
				index_count += 1
				#store the size of the k(:) array
				self.coolZ_nkrates = index_count
				#build reverse as Rji = Rij*gi/gj*exp(-deltaE/T)					
				rate_comment = "!"+str(r_data["up"])+"<-"+str(r_data["down"])+", "+metal_name+" - "+r_data["collider"]+"\n"
				myrate = "k("+str(r_data["rate"])+") * "+str(float(g_up)/float(g_down))+"d0 * exp(-"+deltaE+" * invT)"
				self.coolZ_rates.append(rate_comment+"k("+str(index_count)+") = " + myrate)
				#store the data for the reverse reaction
				rate_data_rev[trans_name_rev] = {"rate": index_count, "collider":r_data["collider"], "up":r_data["down"],\
					"down":r_data["up"]}

			#merge excitation and de-excitation dictionary
			rate_data = dict(rate_data.items() + rate_data_rev.items())
			

			idx_linear_dep_level = 0 #ground level will be removed (for linear dependency)

			#PART 2.2: prepare A matrix for Ax=B system
			#loop on the number of levels (k index, lines of the A matrix)
			nlev = max(level_list)
			collider_list = [] #list of the collider found
			Amatrix = [["0d0" for i in range(nlev+1)] for j in range(nlev+1)] #init matrix to 0d0
			for klev in level_list[:]:
				#if level is suitable for linear depenency use as conservation equation
				if(klev==idx_linear_dep_level):
					Amatrix[klev] = ["1d0" for i in range(nlev+1)]
					continue
				#loop on the number of levels (i index, i<->j)
				for ilev in level_list:
					#loop on the number of levels (j index, i<->j)						
					for jlev in level_list:
						#no transitions from the same level
						if(ilev==jlev): continue
						#k->j transitions (de-populate k level: Akj, Ckj)
						if(klev==ilev):
							A_name = "A("+str(klev+1)+","+str(ilev+1)+")" #Ax=b matrix element name
							trans_name = str(ilev)+"->"+str(jlev) #transition key k->j
							#compute excitation Cij
							#loop on rate data to find keys that contain the transition name (rates k->j)
							for r_key in rate_data:
								r_key_part = r_key.split("_")[0]
								if(trans_name==r_key_part):
									matrix_rate = " &\n- k("+str(rate_data[r_key]["rate"])+") * coll_"+\
										self.convertMetal2F90(rate_data[r_key]["collider"])
									Amatrix[klev][ilev] += matrix_rate
									collider_list.append(rate_data[r_key]["collider"])
							#search transitions k->j (Aij)
							if(trans_name in trans_data):
								Aij_fmt = ("%e" % trans_data[trans_name]["Aij"]).replace("e","d")
								matrix_rate = " &\n- " + Aij_fmt
								Amatrix[klev][ilev] += matrix_rate
								#if photoinduced needed compute it
								if(self.usePhotoInduced and trans_data[trans_name]["Bij"]>0e0):
									Bij_fmt = ("%e" % trans_data[trans_name]["Bij"]).replace("e","d")
									de_eVs = ("%e" % trans_data[trans_name]["denergy_eV"]).replace("e","d")
									matrix_rate = " &\n- "+Bij_fmt+" * get_photoIntensity("+de_eVs+")"
									Amatrix[klev][ilev] += matrix_rate

						#i->k transitions (populate k level: Aik, Cik)
						if(klev==jlev):
							A_name = "A("+str(klev+1)+","+str(ilev+1)+")" #Ax=b matrix element name
							trans_name = str(ilev)+"->"+str(jlev) #transition key i->k
							#add collision, i.e. Cij
							#loop on rate data to find keys that contain the transition name (rates i->k)
							for r_key in rate_data:
								r_key_part = r_key.split("_")[0] #key without _collider (e.g. 2->1)
								if(trans_name==r_key_part):
									matrix_rate = " &\n+ k("+str(rate_data[r_key]["rate"])+") * coll_"+\
										self.convertMetal2F90(rate_data[r_key]["collider"])
									Amatrix[klev][ilev] += matrix_rate
									collider_list.append(rate_data[r_key]["collider"])
							#add transition, i.e. Aij
							#search transitions i->k
							if(trans_name in trans_data):
								Aij_fmt = ("%e" % trans_data[trans_name]["Aij"]).replace("e","d")
								matrix_rate = " &\n+ "+Aij_fmt
								Amatrix[klev][ilev] += matrix_rate
								#if photoinduced needed compute it
								if(self.usePhotoInduced and trans_data[trans_name]["Bij"]>0e0):
									Bij_fmt = ("%e" % trans_data[trans_name]["Bij"]).replace("e","d")
									de_eVs = ("%e" % trans_data[trans_name]["denergy_eV"]).replace("e","d")
									matrix_rate = " &\n+ "+Bij_fmt+" * get_photoIntensity("+de_eVs+")"
									Amatrix[klev][ilev] += matrix_rate

			#PART 2.3: prepare the cooling
			full_B_vector_cool = []
			full_B_vector_heat = []
			for k,t_data in trans_data.iteritems():
				Aij_fmt = ("%e" % t_data["Aij"]).replace("e","d") #f90ish format for Aij
				deltaE = t_data["denergy_K"]
				deltaE_fmt = ("%e" % deltaE).replace("e","d") #f90ish format for deltaE
				photoIB = ""
				if(self.usePhotoInduced and t_data["Bij"]>0e0):
					Bij_fmt = ("%e" % t_data["Bij"]).replace("e","d")
					de_eVs = ("%e" % t_data["denergy_eV"]).replace("e","d")
					photoIB = " &\n + "+Bij_fmt+" * get_photoIntensity("+de_eVs+")"
				#append cooling to final sum over transitions list
				full_B_vector_cool.append("B("+str(t_data["up"]+1)+") * ("+Aij_fmt + photoIB +") * "+deltaE_fmt)

				#add induced heating if needed
				if(self.usePhotoInduced and t_data["Bij"]>0e0):
					Bji_fmt = ("%e" % t_data["Bji"]).replace("e","d")
					de_eVs = ("%e" % t_data["denergy_eV"]).replace("e","d")
					photoIB = Bji_fmt+" * get_photoIntensity("+de_eVs+")"
					full_B_vector_heat.append("B("+str(t_data["up"]+1)+") * "+photoIB +" * "+deltaE_fmt)

			#join the cooling vector as a sum
			full_B_vector = (" &\n + ".join(full_B_vector_cool))
			if(len(full_B_vector_heat)>0):
				full_B_vector += " &\n - " + (" &\n - ".join(full_B_vector_heat))


			#uniqe collider_list
			ucollider_list = []
			for x in collider_list:
				#replace for f90-friendly name
				x = x.replace("+","j").replace("-","k")
				if(x in ucollider_list): continue
				ucollider_list.append(x)
			collider_list = ucollider_list


			#PART 2.4: prepare the function for cooling
			nlev = max(levels_data)+1 #maximum number of levels
			function_name = "cooling"+metal_name #name of the function
			full_function = "!************************************\n"
			full_function += "function "+function_name+"(n,inTgas,k)\n"
			full_function += "use krome_commons\n"
			full_function += "use krome_photo\n"
			full_function += "implicit none\n"

			#declaration of varaibles
			full_function += "integer::i, hasnegative, nmax\n"
			full_function += "real*8::"+function_name+",n(:),inTgas,k(:)\n"
			full_function += "real*8::A("+str(nlev)+","+str(nlev)+"),Ain("+str(nlev)+","+str(nlev)+")\n"
			full_function += "real*8::B("+str(nlev)+"),tmp("+str(nlev)+")\n"

			#declaration of alias variable for colliders
			for x in collider_list:
				full_function += "real*8::coll_"+x+"\n"

			full_function += "\n!colliders should be >0\n"

			#initialization of colliders (replace ortho/para H2)
			for x in collider_list:
				xn = "n(idx_"+x+")"
				if(x=="H2or"): xn = "n(idx_H2) * phys_orthoParaRatio / (phys_orthoParaRatio+1d0)"
				if(x=="H2pa"): xn = "n(idx_H2) / (phys_orthoParaRatio+1d0)"
				full_function += "coll_"+x+" = max("+xn+", 0d0)\n" #collider must be positive

			full_function += "\n!deafault cooling value\n"
			full_function += function_name +" = 0d0\n\n" #default
			full_function += "if(n(idx_"+metal_name_f90+")<1d-15) return\n\n" #if low coolant abundance skip all
			full_function += "A(:,:) = 0d0\n\n" #init A matrix to zero

			#write the initialization of first column of the A matrix 
			# (will be used by the f90 to reduce the size of the problem)
			for j in range(nlev):
				if(Amatrix[j][0]!="0d0"):
					matrix_element = Amatrix[j][0].replace("0d0 &\n","")
					full_function += "A("+str(j+1)+",1) = "+matrix_element+"\n"

			#write the initialization of diagonal elements of the A matrix 
			# (will be used by the f90 to reduce the size of the problem)
			for j in range(1,nlev):
				if(Amatrix[j][j]!="0d0"):
					matrix_element = Amatrix[j][j].replace("0d0 &\n","")
					full_function += "A("+str(j+1)+","+str(j+1)+") = "+matrix_element+"\n"


			#the size of the problem can be reduced up to the last non-zero row of the left triangular matrix
			# only interesting with more levels
			if(nlev>3):
				full_function += "\n!reduce the size of the problem if possible\n" #reverse is faster
				full_function += "nmax = 1\n"
				full_function += "do i="+str(nlev)+",2,-1\n"
				full_function += " if(A(i,1)>0d0) then\n"
				full_function += "  nmax = i\n" #store new size of the problem
				full_function += "  exit\n" #break loop
				full_function += " end if\n"
				full_function += "end do\n\n"

				#control if the problem is not 1-level
				full_function += "!no need to solve a 1-level problem\n"
				full_function += "if(nmax==1) return\n\n"


			#write the A matrix column-wise. A(:,1) matrix column computed above
			# as well as the diagonal elements 
			for i in range(1,nlev):
				for j in range(nlev):
					#skip diagonal elements since already written (see above)
					if(i==j): continue
					#skip zero elements
					if(Amatrix[j][i]!="0d0"):
						matrix_element = Amatrix[j][i].replace("0d0 &\n","")
						full_function += "A("+str(j+1)+","+str(i+1)+") = "+matrix_element+"\n"

			full_function += "\n!build matrix B\n"
			full_function += "B(:) = 0d0\n" #default B matrix
			full_function += "B("+str(idx_linear_dep_level+1)+") = n(idx_"+metal_name_f90+")\n\n" #conservation equation RHS
			full_function += "Ain(:,:) = A(:,:)\n\n" #store initial matrix for debug purposes


			#choose the correct solver depending on the number of levels
			if(nlev>3):
				full_function += "call mydgesv(nmax, A(:,:), B(:), \""+function_name+"\")\n\n"
				self.needLAPACK = True
			elif(nlev==2):
				full_function += "call mylin2(A(:,:), B(:))\n\n"
			elif(nlev==3):
				full_function += "call mylin3(A(:,:), B(:))\n\n"
			else:
				sys.exit("ERROR: strange number of levels for linear system in Zcooling: "+str(nlev))

			full_function += "!store population\n"
			full_function += "pop_level_"+metal_name+"(:) = B(:)\n"

			#negative small values can be flushed to 1d-40
			full_function += "!sanitize negative values\n"
			full_function += "hasnegative = 0\n"
			full_function += "do i=1,"+str(nlev)+"\n"
			full_function += " if(B(i)<0d0) then\n"
			full_function += "  if(abs(B(i)/n(idx_"+metal_name_f90+"))>1d-10) then\n"
			full_function += "   hasnegative = 1\n"
			full_function += "  else\n"
			full_function += "   B(i) = 1d-40\n"
			full_function += "  end if\n"
			full_function += " end if\n"
			full_function += "end do\n\n"

			#when negative value are found print some stuff and stop
			full_function += "!check if B has negative values\n"
			full_function += "if(hasnegative>0)then\n"
			full_function += " print *,\"ERROR: minval(B)<0d0 in "+function_name+"\"\n"
			full_function += " print *,\"ntot_"+metal_name+" =\", n(idx_"+metal_name_f90+")\n"
			full_function += " print *,\"Tgas =\", inTgas\n"
			full_function += " print *,\"B(:) unrolled:\"\n"
			full_function += " do i=1,size(B)\n"
			full_function += "  print *, i, B(i)\n"
			full_function += " end do\n"
			full_function += " print *,\"A(:,:) min/max:\"\n"
			full_function += " do i=1,size(B)\n"
			full_function += "  print *, i, minval(Ain(i,:)), maxval(Ain(i,:))\n"
			full_function += " end do\n\n"
			full_function += " print *,\"A(:,:)\"\n"
			full_function += " do i=1,size(B)\n"
			full_function += "  tmp(:) = Ain(i,:)\n"
			full_function += "  print '(I5,99E17.8)', i, tmp(:)\n"
			full_function += " end do\n"
			full_function += " stop\n"
			full_function += "end if\n\n"

			#when the population for each level is known compute the cooling (see above) 
			full_function += function_name + " = " +full_B_vector+"\n\n"
			full_function += "end function "+function_name+"\n\n"

			#append the function to the list of the functions
			self.coolZ_functions.append([function_name,full_function])
			self.coolZ_poplevelvars.append("pop_level_"+metal_name+"("+str(nlev)+")")

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
					fout.write("\tinteger,parameter::nPhotoBins=" + str(self.photoBins) + "\n")
					fout.write("\tinteger,parameter::nPhotoRea=" + str(self.nPhotoRea) + "\n")

			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			elif(srow == "#KROME_phys_commons"):
				for x in self.physVariables:
					fout.write("real*8::phys_"+x[0]+"\n")
			elif(srow == "#KROME_cool_index"):
				idxcool = get_cooling_index_list()
				for x in idxcool:
					fout.write("integer,parameter::"+x+"\n")
			elif(srow == "#KROME_heat_index"):
				idxheat = get_heating_index_list()
				for x in idxheat:
					fout.write("integer,parameter::"+x+"\n")
			elif(srow == "#KROME_implicit_arr_r"):
				for j in range(self.maxnreag):
					fout.write("integer::arr_r"+str(j+1)+"(nrea)\n")
			elif(srow == "#KROME_implicit_arr_p"):
				for j in range(self.maxnprod):
					fout.write("integer::arr_p"+str(j+1)+"(nrea)\n")
			elif(srow == "#KROME_photoheating_variables" and self.useHeatingPhoto):
				fout.write("real*8::"+(",".join(pheatvars))+"\n")
			elif(srow == "#KROME_opt_variables"):
				fout.write(optVariables)
			elif(srow == "#KROME_user_commons"):
				if(len(self.commonvars)>0):
					fout.write("real*8::"+(",".join(self.commonvars))+"\n")
					fout.write("!$omp threadprivate("+(",".join(self.commonvars))+")\n")
			elif(srow == "#KROME_photobins_array"):
				if(self.photoBins>0):
					fout.write("real*8::photoBinJ(nPhotoBins) !intensity per bin, erg/s/sr/Hz/cm2\n")
					fout.write("real*8::photoBinEleft(nPhotoBins) !left limit of the freq bin, eV\n")
					fout.write("real*8::photoBinEright(nPhotoBins) !right limit of the freq bin, eV\n")
					fout.write("real*8::photoBinEmid(nPhotoBins) !middle point of the freq bin, eV\n")
					fout.write("real*8::photoBinEdelta(nPhotoBins) !size of the freq bin, eV\n")
					fout.write("real*8::photoBinEidelta(nPhotoBins) !inverse of the size of the freq bin, 1/eV\n")
					fout.write("real*8::photoBinJTab(nPhotoRea,nPhotoBins) !xsecs table, cm2\n")
					fout.write("real*8::photoBinRates(nPhotoRea) !photo rates, 1/s\n")
					fout.write("real*8::photoBinHeats(nPhotoRea) !photo heating, erg/s\n")
					fout.write("real*8::photoBinEth(nPhotoRea) !energy treshold, eV\n")
					fout.write("!$omp threadprivate(photoBinJ,photoBinEleft,photoBinEright,photoBinEmid,photoBinEdelta, &\n")
					fout.write("!$omp               photoBinEidelta,photoBinJTab,photoBinRates,photoBinHeats,photoBinEth)\n")
			#write the anytab common variables
			elif(srow == "#KROME_vars_anytab"):
				stab = ""
				for i in range(len(self.anytabvars)):
					tabvar = self.anytabvars[i]
					tabfile = self.anytabfiles[i]
					tabsize = self.anytabsizes[i]
					stab += "real*8::" + tabvar+"_anytabx("+tabsize[0]+")\n"
					stab += "real*8::" + tabvar+"_anytaby("+tabsize[1]+")\n"
					stab += "real*8::" + tabvar+"_anytabz("+tabsize[0]+","+tabsize[1]+")\n"
					stab += "real*8::" + tabvar+"_anytabxmul\n"
					stab += "real*8::" + tabvar+"_anytabymul\n"
				fout.write(stab+"\n")

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
		constants.append(["iplanck_eV","1d0/planck_eV","1 / eV / s"])
		constants.append(["iplanck_J","1d0/planck_J","1 / J / s"])
		constants.append(["iplanck_erg","1d0/planck_erg","1 / erg / s"])
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
		constants.append(["hubble","0.704d0","dimensionless"])
		constants.append(["Omega0","1.0d0","dimensionless"])
		constants.append(["Omegab","0.0456d0","dimensionless"])
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

		#preapare metallicity functions
		#metallicity dictionary
		metalDict = dict()
		for x in specs:
			for atom,count in x.atomcount2.iteritems():
				if(atom in metalDict):
					metalDict[atom].append([x.fidx,count])
				else:
					metalDict[atom] = [[x.fidx,count]]

		#metallicity functions
		nZs = ["H","He","+","-","E"] #these are not metals
		zGets = []
		for k,v in metalDict.iteritems():
			if(k in nZs): continue #skip non-metals
			parts = []
			for x in v:
				smult = "" #multiplication factor string if >1 atom/particle
				if(x[1]>1): smult = str(x[1])+"d0*"
				parts.append(smult+"n("+x[0]+")")
			zGets.append([k, "z"+k, (" &\n + ".join(parts))])
				

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
                        if(srow == "#IFKROME_useXrays" and not(self.useXRay)): skip = True

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
			#replace arrays for best flux
			elif(srow == "#KROME_arr_reactprod"):
				for i in range(self.maxnreag):
					fout.write("if(arr_r"+str(i+1)+"(i) == idx_found) found = .true.\n") 
				for i in range(self.maxnprod):
					fout.write("if(arr_p"+str(i+1)+"(i) == idx_found) found = .true.\n") 
			elif(srow == "#KROME_conserve"):
				fout.write(krome_conserve+"\n")
			elif(srow == "#KROME_col2num_method"):
				if(self.columnDensityMethod=="DEFAULT"):
					fout.write("col2num = 1d3 * (ncalc/1.8d21)**1.5\n")
				elif(self.columnDensityMethod=="JEANS"):
					fout.write("col2num = 2d0 * ncalc / get_jeans_length(n(:),Tgas)\n")
				else:
					sys.exit("ERROR: method "+self.columnDensityMethod+" unknown for col2num")
			elif(srow == "#KROME_num2col_method"):
				if(self.columnDensityMethod=="DEFAULT"):
					fout.write("num2col = 1.8d21*(max(ncalc,1d-40)*1d-3)**(2./3.)\n")
				elif(self.columnDensityMethod=="JEANS"):
					fout.write("num2col = 0.5d0 * ncalc * get_jeans_length(n(:),Tgas)\n")
				else:
					sys.exit("ERROR: method "+self.columnDensityMethod+" unknown for num2col")
			elif(srow == "#KROME_metallicity_functions"):
				solar = get_solar_abundances() #get solar abundances
				ffs = "" #metallicity functions
				for zg in zGets:
					if(not(zg[0] in solar)): continue #skip if solar abundance is not present
					#prepare function
					ffname = "get_metallicity"+zg[0]
					ff = "!*****************************\n"
					ff += "! get metallicity using "+zg[0]+" as reference\n"
					ff += "function "+ffname+"(n)\n"
					ff += "use krome_commons\n"
					ff += "implicit none\n"
					ff += "real*8::n(:),"+ffname+","+zg[1]+",nH\n"
					ff += "nH = get_Hnuclei(n(:))\n"
					ff += zg[1]+" = "+zg[2]+"\n"
					ff += zg[1]+" = max("+zg[1]+", 0d0)\n\n"
					ff += ffname + " = log10("+zg[1]+"/nH+1d-40) - "+solar[zg[0]]+"\n\n" #compute metallicity
					ff += "end function "+ffname+"\n\n"
					ffs += ff #append function to the others
				fout.write(ffs+"\n") #replace functions
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
				#look for the largest qeff value
				maxqeff = 0e0
				for x in reacts:
					maxqeff = max(maxqeff,x.qeff)
				#if 0e0 is the largest then compress the array, else write it explicitely
				if(maxqeff==0e0):
					fout.write("get_qeff(:) = 0e0\n")
				else:
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
					if(not("H" in x.atomcount2)): continue
					if(x.atomcount2["H"]==0): continue
					hmult = ("*"+format_double(x.atomcount2["H"]) if x.atomcount2["H"]>1 else "")
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

				elif(self.typeGamma=="EXACT" or self.typeGamma=="VIB" or self.typeGamma=="ROT" or self.typeGamma=="REDUCED"):
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
							#skip every diatoms except H2 and CO if REDUCED
							if(gtype=="REDUCED" and (mol.name!="H2" and mol.name!="CO")): continue 
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

		#replace photoionization and photoheating functions
		skip = skip_heat = False
		for row in fh:
			srow = row.strip()
			if(row.strip() == "#IFKROME_usePhIoniz" and not(self.usePhIoniz)): skip = True
			if(row.strip() == "#IFKROME_usePhotoBins" and not(self.photoBins>0)): skip = True
			if(row.strip() == "#ENDIFKROME"): skip = False

			if(row.strip() == "#IFKROME_photobin_heat" and not(self.useHeatingPhoto)): skip_heat = True
			if(row.strip() == "#ENDIFKROME_photobin_heat"): skip_heat = False

			if(skip or skip_heat): continue

			#replace pragma with the initialization of the photorate table in bins
			if(srow=="#KROME_photobin_xsecs"):
				phbinx = ""
				for rea in reacts:
					if(rea.kphrate==None): continue
					phbinx += "\n!"+rea.verbatim+"\n"
					phbinx += "kk = "+rea.kphrate+"\n"
					phbinx += "if(energy_eV<"+str(rea.Tmin)+") kk = 0d0\n"
					phbinx += "if(energy_eV>"+str(rea.Tmax)+") kk = 0d0\n"
					phbinx += "photoBinJTab("+str(rea.idxph)+",j) = kk\n"
				row = phbinx+"\n"
			#replace the energy treshold assuming that it is equal to Tmin
			elif(srow=="#KROME_photobin_Eth"):
				phbinx = ""
				for rea in reacts:
					if(rea.kphrate==None): continue
					phbinx += "photoBinEth("+str(rea.idxph)+") = "+str(rea.Tmin)+" !"+rea.verbatim+"\n"
				row = phbinx+"\n"
			#replace pragma with the opacity calculation as N_i*sigma_i for any species
			elif(srow=="#KROME_photobin_opacity" and self.usePhotoOpacity):
				phbintau = ""
				#loop on the species looking for photorates
				for rea in reacts:
					if(rea.kphrate==None): continue
					phbintau += "tau = tau + photoBinJTab("+str(rea.idxph)+",j) * ncol("+rea.reactants[0].fidx+") !"\
						+rea.verbatim+"\n"
				row = phbintau+"\n"

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
		if(countNoTab>0 and self.useTabs):
			if(len(coevars)!=0):
				#define variables
				kvars = "real*8::"+(",".join([x.strip() for x in coevars.keys()]))
				#variables initialization
				klist = [[k+" = "+v[1]+"\n",v[0]] for k,v in coevars.iteritems()] #this mess is to sort dict
				klist = sorted(klist, key=lambda x: x[1])
				klist = "".join([x[0] for x in klist])

		#prepares the reaction modifiers
		#tokenize to replace k(:) with coe_tab(:)
		import tokenize,cStringIO
		kModifierFull = "" #string that will contans all the lines of the rate modifier
		for kmod in self.kModifier:
			#string to tokenizable object
			src = cStringIO.StringIO(kmod).readline
			#tokenize
			tokenized = tokenize.generate_tokens(src)
			kmodTok = "" #string with the k->coe_tab replaced
			for tok in tokenized:
				if(tok[1]=="k"):
					kmodTok += "coe_tab"
				else:
					kmodTok += tok[1]
			#comments are not tokenized
			if(kmod[0]=="!"): kmodTok = kmod
			#append the correct string
			kModifierFull += kmodTok+"\n"


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
			row = row.replace("#KROME_rateModifier", kModifierFull)

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
		from math import sqrt
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


		#bremsstrahlung (free-free) for all the ions as charge**2*n_ion
		bms_ions = "bms_ions ="
		#look for ions (charge>0)
		for x in specs:
			charge = x.charge #store species charge
			if(not(x.is_atom)): continue #only atoms for free-free
			if(charge>0):
				mult = "" #multiplication factor
				if(charge>1): mult = str(charge*charge)+".d0*" #multipy by the square of the charge
				bms_ions += " +"+mult+"n("+x.fidx+")" #add Z^2*abundance


		#load data for free-bound
		fhfb = open("data/ip.dat","rb")
		fbdata = []
		for row in fhfb:
			srow = row.strip()
			if(srow==""): continue
			if(srow[0]=="#"): continue
			arow = [x for x in srow.split(" ") if x!=""]
			#Z: atomic number, ion: ionization degree (e.g. HII=1), energy_eV: ioniz potential, n0: principal quantum number 
			fbdata.append({"Z":int(arow[0]), "ion":int(arow[1]), "energy_eV":float(arow[5]), "n0":int(arow[6])})

		skip = skip_nleq = False
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
			if(srow == "#IFKROME_useCoolingFF" and not(self.useCoolingFF)): skip = True
			if(srow == "#IFKROME_useCoolingContinuum" and not(self.useCoolingCont)): skip = True
			if(srow == "#IFKROME_useLAPACK" and not(self.needLAPACK)): skip = True #skip calls to LAPACK
			if(srow == "#IFKROME_use_NLEQ" and not(self.useNLEQ)): skip_nleq = True #skip calls to NLEQ

			if(srow == "#ENDIFKROME_use_NLEQ"): skip_nleq = False
			if(srow == "#ENDIFKROME"): skip = False

			if(skip or skip_nleq): continue

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
			elif(row.strip() == "#KROME_escape_vars"):
				if(len(self.fcn_levs)>0):
					fcn_levs = sorted(self.fcn_levs)
					bvars = "real*8::"+(",".join(["beta"+str(x)+"("+str(x)+","+str(x)+")" for x in fcn_levs]))+"\n"
					bvars += "real*8::"+(",".join(["Besc"+str(x)+"("+str(x)+","+str(x)+")" for x in fcn_levs]))+"\n"
					bvars += "real*8::"+(",".join(["MMesc"+str(x)+"("+str(x)+","+str(x)+")" for x in fcn_levs]))+"\n"
					bvars += "real*8::ntotcoll\n"
					fout.write(bvars)
			elif(row.strip() == "#KROME_fcn_cases"):
				if(len(self.fcn_levs)>0):
					fcn_levs = sorted(self.fcn_levs)
					for i in range(len(fcn_levs)):
						preif = ""
						if(i>0): preif = "else "
						fcase = preif+"if(n=="+str(fcn_levs[i])+") then\n"		
						fcase += "call fcn_"+str(fcn_levs[i])+"(n,x(:),f(:))\n"
						fout.write(fcase)
					fout.write("else\n")
					fout.write("print *,\"ERROR: unknown case in fcn subroutine!\", n\n")
					fout.write("stop\n")
					fout.write("end if\n")
			#prepare free-bound calculation
			elif(row.strip()=="#KROME_FB_cooling_data"):
				for x in self.specs:
					charge = x.charge #species charge
					Zatom = x.zatom #species atomic number
					if(not(x.is_atom)): continue #skip molecules
					if(charge<=0): continue #skip neutral and anions
					#search the data in the table loaded above
					dataFound = False
					for dataip in fbdata:
						if(dataip["Z"]==Zatom and dataip["ion"]==charge):
							mydataip = dataip #store the line
							dataFound = True #update found flag
							break #break the loop 
					#raise error if no data found
					if(not(dataFound)): sys.exit("ERROR: no data found for "+x.name)
					#start calculation (precompute most of the known stuff)
					n0 = mydataip["n0"] #principal quantum number
					ion = mydataip["ion"] #ionization degree (charge)
					E_eV = mydataip["energy_eV"] #ionization potential eV
					z0 = sqrt(E_eV/1.359808e1)*n0
					E0 = E_eV/1e3 #keV
					En1 = (ion)**2*1.359808e1/(n0+1)**2/1e3 #keV
					l0 = 911.9e0*(n0/z0)**2 #angstrom
					ln1 = 911.9e0/ion**2*(float(n0+1))**2 #angstrom
					dzion = Zatom-ion
					if(dzion>22):
						zeta0 = -dzion+55e0
					elif(dzion<=22 and dzion>8):
						zeta0 = -dzion+27e0
					elif(dzion<=8 and dzion>0):
						zeta0 = -dzion+9e0
					else:
						zeta0 = -dzion+1e0
					print x.name
					print "f2 = "+format_double(.9*zeta0*z0**4/n0**5)+"*exp("+format_double(0.1578e0*z0**2/n0**2)+"*invT6)"
					print "gfb = "+"0.1578d0*n("+x.fidx+")*f2*invT6"

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
					fout.write(x+"\n\n")
			elif(row.strip() == "#KROME_coolingZ_popvars"):
				if(len(self.coolZ_poplevelvars)>0):
					for popvar in self.coolZ_poplevelvars:
						fout.write("real*8::"+popvar+"\n")
					for popvar in self.coolZ_poplevelvars:
						funct_name = popvar.split("(")[0]
						fout.write("!$omp threadprivate("+funct_name+")\n")
			elif(row.strip() == "#KROME_popvar_dump"):
				if(len(self.coolZ_poplevelvars)>0):
					for popvar in self.coolZ_poplevelvars:
						funct_name = popvar.split("(")[0]
						metal_name = funct_name.split("_")[-1]
						fout.write("!"+popvar+"\n")
						fout.write("do i=1,size("+funct_name+")\n")
						fout.write(" write(nfile,'(a8,I5,2E17.8e3)') \""+metal_name+"\", i, Tgas, "+funct_name+"(i)\n")
						fout.write("end do\n\n")
			else:
				#replace pragma for total metals
				row = row.replace("#KROME_tot_metals", self.totMetals)
				
				if(self.H2opacity=="RIPAMONTI"):
					#thick case (note that 1.25d-10 = 1/8e9)
					row = row.replace("#KROME_H2opacity", "&\n* min(1.d0, max(1.25d-10 * sum(n(1:nmols)),1d-40)**(-.45))")
				elif(self.H2opacity=="OMUKAI"):
					#thick case using table provided by Omukai (priv. comm. 2014)
					row = row.replace("#KROME_H2opacity", "&\n* H2opacity_omukai(Tgas, n(:))")
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
				#phstuff = get_ph_stuff(react)
				#if(phstuff==None): continue
				#reaname = phstuff["reaname"]
				#reag = react.reactants
				#fake_opacity = ""
				#if(self.useFakeOpacity): fake_opacity = " * exp(-n(" + reag[0].fidx + ") / n0)"
				if(react.idxph<=0): continue
				pheatvars.append("photoBinHeats("+str(react.idxph)+") * n(" + react.reactants[0].fidx + ")")

		#replace pragma with strings built above
		skip = False
		for row in fh:
			srow = row.strip()
			if(row.strip() == "#KROME_header"):

				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			else:
				if(row.strip() == "#IFKROME_useHeatingCR" and not(self.useHeatingCR)): skip = True
				if(row.strip() == "#IFKROME_useHeatingdH" and (not(self.useHeatingdH) or len(dH_varsa)==0)): skip = True
				if(row.strip() == "#IFKROME_useHeatingCompress" and not(self.useHeatingCompress)): skip = True
				if(row.strip() == "#IFKROME_useHeatingPhoto" and not(self.useHeatingPhoto)): skip = True
				if(row.strip() == "#IFKROME_useHeatingPhotoAv" and not(self.useHeatingPhotoAv)): skip = True
				if(row.strip() == "#IFKROME_useHeatingPhotoDust" and not(self.useHeatingPhotoDust)): skip = True
				if(row.strip() == "#IFKROME_useHeatingXRay" and not(self.useHeatingXRay)): skip = True
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

				#replace cosmic ray heating
				if("#KROME_heatingCR" in srow):
					CRheat = ""
					for rea in self.reacts:
						if(not(rea.isCR)): continue
						CRheat += "!"+rea.verbatim+"\n"
						CRheat += "heat_CR = heat_CR + k("+str(rea.idx)+") * n("+rea.reactants[0].fidx+") * Hfact\n\n"
					row = row.replace("#KROME_heatingCR",CRheat)

				if(len(pheatvars)>0):
					row = row.replace("#KROME_photo_heating", "photo_heating = " + (" &\n+ ".join(pheatvars)))
				row = row.replace("#KROME_HChem_terms", HChem) #replace chemical heating terms
				row = row.replace("#KROME_HChem_dust", HChemDust) #replace chemical heating for dust

				#replace metallicity
				if("#KROME_photoDustZ" in row):
					zFound = False
					for zz in ["Fe","C","O","Si"]:
						for x in self.specs:
							if(zz==x.name):
								zFound = True
								dustZ = zz
								break
						if(zFound): break
					if(not(zFound)): 
						row = row.replace("#KROME_photoDustZ","0d0")
					else:
						row = row.replace("#KROME_photoDustZ","1d1**get_metallicity"+zz+"(n(:))")

				#replace correct dissociation rates
				if("#KROME_RdissH2" in row):
					rdh2Found = False
					for rea in self.reacts:
						R = sorted([x.name for x in rea.reactants])
						P = sorted([x.name for x in rea.products])
						if(R==["H2"] and P==["H","H"]):
							rateDissH2 = "k("+str(rea.idx)+")"
							rdh2Found = True
							break
					#check if rate photodissiocation rate is present in the network
					if(not(rdh2Found)): 
						print "ERROR: if you use PHOTOAV heating you should have"
						print " H2 photodissiocation rate in your chemical network!"
						sys.exit()			

					row = row.replace("#KROME_RdissH2",rateDissH2) #replace pragma with H2 photodissociation rate

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
			if(srow == "#IFKROME_use_thermo_toggle" and not(self.useThermoToggle)): skip = True
			if(srow == "#IFKROME_report" and not(self.doReport)): skip = True
			if(srow == "#IFKROME_useDust" and not(self.useDust)): skip = True

			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			coolPragmaFound = False
			#include cooling cmb floor if necessary
			if("#KROME_coolCMBfloor" in srow):
				coolPragmaFound = True
				if(self.useCoolCMBFloor):
					srow = srow.replace("#KROME_coolCMBfloor"," + cooling(n(:), phys_Tcmb)")
				else:
					srow = srow.replace("#KROME_coolCMBfloor","")
				
			#replace quenching function for cooling
			if("#KROME_coolingQuench" in srow):
				coolPragmaFound = True
				if(self.coolingQuench<0e0):
					srow = srow.replace("#KROME_coolingQuench","")
				else:
					qfunc = " &\n * 0.5d0 * (tanh(Tgas - "+format_double(self.coolingQuench)+") + 1d0)"
					srow = srow.replace("#KROME_coolingQuench",qfunc)

			#quench and cmbfloor are on the same line so write the replacements togheter
			if(coolPragmaFound):
				fout.write(srow+"\n")
				continue

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

			#replace the pragma with the computation of the photorates using the opacity computed with
			# the approximation of Glover+2009 Eqn.2
			elif(srow == "#KROME_photobins_compute_thick" and self.usePhotoOpacity):
				fout.write("call calc_photoBins_thick(n(:))\n")
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

			elif(srow == "#KROME_ODEModifier"):
				#write the ODE modifiers
				odeModifierFull = "" #string that will contans all the lines of the ode modifier
				for kmod in self.odeModifier:
					#append the correct string
					fout.write(kmod+"\n")

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

		solar = get_solar_abundances()

		scaleZ = []
		#looks for H to rescale the metallicity otherwise skips
		has_H = False
		sHtot = "Htot = "
		for mols in specs:
			if(not("H" in mols.atomcount2)): continue 
			Hcount = mols.atomcount2["H"]
			if(Hcount>0): has_H = True
			if(Hcount==1): sHtot += " &\n + n("+mols.fidx+")"
			if(Hcount>1): sHtot += " &\n + "+str(Hcount)+"d0 * n("+mols.fidx+")"

		scaleZ.append(sHtot) #Htot= is the first of the list 
		#creates the metallicity rescaling subroutine
		for (k,v) in solar.iteritems():
			if(not(has_H)):
				scaleZ = [] #reset scaleZ since Htot= is no longer needed
				break #skip routine if H is not present
			for mols in specs:
				if(k.upper()==mols.name.upper()):
					scaleZ.append("n("+mols.fidx+") = max(Htot * 1d1**(Z+log10("+str(v)+")), 1d-40)")

		for row in fh:

			srow = row.strip()

			if(srow == "#IFKROME_usePhotoBins" and self.photoBins<=0): skip = True
			if(srow == "#IFKROME_useStars" and not(self.useStars)): skip = True
			if(srow == "#IFKROME_use_cooling" and not(self.use_cooling)): skip = True
			if(srow == "#IFKROME_use_thermo" and not(self.use_thermo)): skip = True
			if(srow == "#IFKROME_use_coolingZ" and not(self.useCoolingZ)): skip = True
			if(srow == "#IFKROME_useXrays" and not(self.useXRay)): skip = True

			if(srow == "#ENDIFKROME"): skip = False

			if(skip): continue

			if(srow == "#KROME_species"):
				for x in specs:
					fout.write("\tinteger,parameter::" + "KROME_"+x.fidx + " = " + str(x.idx) +"\t!"+x.name+"\n")
			elif(srow == "#KROME_user_commons_functions"):
				funcs = ""
				for x in self.commonvars:
					fsetname = "krome_set_"+x
					fset = "\n!*******************\n"
					fset += "subroutine "+fsetname+"(argset)\n"
					fset += "use krome_commons\n"
					fset += "implicit none\n"
					fset += "real*8::argset\n"
					fset += x+" = argset\n"
					fset += "end subroutine "+fsetname+"\n"
					
					fgetname = "krome_get_"+x
					fget = "\n!*******************\n"
					fget += "function "+fgetname+"()\n"
					fget += "use krome_commons\n"
					fget += "implicit none\n"
					fget += "real*8::"+fgetname+"\n"
					fget += fgetname+" = "+x+"\n"
					fget += "end function "+fgetname+"\n"
					funcs += fset + fget
				fout.write(funcs)

			elif(srow == "#KROME_cool_index"):
				idxcool = get_cooling_index_list()
				for x in idxcool:
					fout.write("real*8,parameter::krome_"+x+"\n")

			elif(srow == "#KROME_heat_index"):
				idxheat = get_heating_index_list()
				for x in idxheat:
					fout.write("real*8,parameter::krome_"+x+"\n")
			elif(srow=="#KROME_print_phys_variables"):
					for x in self.physVariables:
						fout.write("print *, \""+x[0]+":\", phys_"+x[0]+"\n")
			elif(srow=="#KROME_set_get_phys_functions"):
				for x in self.physVariables:
					#set subroutine
					funcname = "krome_set_"+x[0]
					fout.write("!*******************\n")
					fout.write("subroutine "+funcname+"(arg)\n")
					fout.write(" use krome_commons\n")
					fout.write(" implicit none\n")
					fout.write(" real*8::arg\n")
					fout.write(" phys_"+x[0]+" = arg\n")
					fout.write("end subroutine "+funcname+"\n\n")

					#get function
					funcname = "krome_get_"+x[0]
					fout.write("!*******************\n")
					fout.write("function "+funcname+"()\n")
					fout.write(" use krome_commons\n")
					fout.write(" implicit none\n")
					fout.write(" real*8::"+funcname+"\n")
					fout.write(funcname+" = phys_"+x[0]+"\n")
					fout.write("end function "+funcname+"\n\n")
			#write the user alias for the cooling functions
			elif(srow == "#KROME_cooling_functions"):
				for x in self.coolZ_functions:
					funcname =  "krome_"+x[0];
					fout.write("\n!*******************\n")
					fout.write("function "+funcname+"(xin,inTgas)\n")
					fout.write("use krome_commons\n")
					fout.write("use krome_cooling\n")
					fout.write("use krome_constants\n")
					fout.write("real*8::xin(:),n(nspec),inTgas,k(nZrate),"+funcname+"\n")
					fout.write("n(:) = 0d0\n")
					fout.write("n(idx_Tgas) = inTgas\n")
					fout.write("n(1:nmols) = xin(:)\n")
					fout.write("k(:) = coolingZ_rate_tabs(inTgas)\n")
					fout.write(funcname+" = "+x[0]+"(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg\n")
					fout.write("end function "+funcname+"\n")

			elif(srow == "#KROME_header"):
				fout.write(get_licence_header(self.version, self.codename,self.shortHead))
			elif(srow == "#KROME_zero_electrons"):
				#check if electron exists
				for x in specs:
					if(x.name=="E"):
						fout.write("x(idx_e) = 0.d0\n")
						break
			elif(srow == "#KROME_electrons_balance"):
				#check if electron exists
				for x in specs:
					if(x.name=="E"):
						fout.write("x(idx_e) = ee\n")
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
				fout.write("\tinteger,parameter::krome_nrea=" + str(self.nrea) + "\n")
				fout.write("\tinteger,parameter::krome_nmols=" + str(nmols) + "\n")
				fout.write("\tinteger,parameter::krome_nspec=" + str(len(specs)) + "\n")
				fout.write("\tinteger,parameter::krome_ndust=" + str(dustArraySize*dustTypesSize) + "\n")
				fout.write("\tinteger,parameter::krome_ndustTypes=" + str(dustTypesSize) + "\n")
				fout.write("\tinteger,parameter::krome_nPhotoBins=" + str(self.photoBins) + "\n")
				fout.write("\tinteger,parameter::krome_nPhotoRates=" + str(self.nPhotoRea) + "\n")
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

			if(srow == "#IFKROME_usePhotoBins" and not(self.photoBins>0)): skip = True
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

			if(srow == "#IFKROME_useH2esc_omukai" and (self.H2opacity!="OMUKAI")): skip = True
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
			#write the anytab initializations
			elif(srow == "#KROME_init_anytab"):
				stab = ""
				for i in range(len(self.anytabvars)):
					tabvar = self.anytabvars[i]
					tabfile = self.anytabfiles[i]
					anytabx = tabvar+"_anytabx(:)"
					anytaby = tabvar+"_anytaby(:)"
					anytabz = tabvar+"_anytabz(:,:)"
					anytabxmul = tabvar+"_anytabxmul"
					anytabymul = tabvar+"_anytabymul"
					stab += "call init_anytab2D(\""+tabfile+"\","+anytabx+",&\n"\
						+anytaby+","+anytabz+",&\n"+anytabxmul+","\
						+anytabymul+")\n"
					stab += "call test_anytab2D(\""+tabvar+"\","+anytabx+",&\n"\
						+anytaby+","+anytabz+",&\n"+anytabxmul+","\
						+anytabymul+")\n"
				fout.write(stab+"\n")

			elif(srow == "#KROME_init_phys_variables"):
				for x in self.physVariables:
					fout.write("phys_"+x[0]+" = "+x[1]+"\n")					
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
		if(self.useDust):
			print "- copying optical data for dust..."
			shutil.copyfile("data/optC.dat", buildFolder+"optC.dat")
			shutil.copyfile("data/optSi.dat", buildFolder+"optSi.dat")

		#copy OMUKAI datafile
		if(self.H2opacity=="OMUKAI"):
			shutil.copyfile("data/escape_H2.dat", buildFolder+"escape_H2.dat")

		#copy file that contains table as indicated by the anytab reactions
		print "- copying anytab files..."
		for i in range(len(self.anytabvars)):
			shutil.copyfile(self.anytabpaths[i], buildFolder+self.anytabfiles[i])

		#copy static files to build
		if(self.is_test):
			print "- copying test to /build..."
			mypath = "tests/"+test_name
			files = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
			for fdir in files:
				shutil.copyfile("tests/"+test_name+"/"+fdir, buildFolder+"/"+fdir)
				print "- copying "+fdir+" to "+buildFolder
			#if Makefile is not present in the tests directory use the default Makefile
			if(not(os.path.exists("tests/"+test_name+"/Makefile"))):
				if(self.useDvodeF90):
					shutil.copyfile("tests/MakefileF90", buildFolder+"Makefile")
				elif(self.buildCompact):
					shutil.copyfile("tests/MakefileCompact", buildFolder+"Makefile")
				else:
					if(self.pedanticMakefile):
						shutil.copyfile("tests/Makefile_pedantic", buildFolder+"Makefile")
					else:
						shutil.copyfile("tests/Makefile", buildFolder+"Makefile")

			#test_file = "tests/"+test_name+"/test.f90"
			#plot_file = "tests/"+test_name+"/plot.gps"

			#chech if test file exists
			#if(not(os.path.isfile(filename))): die("ERROR: Test file \""+test_file+"\" doesn't exist!")

			#shutil.copyfile(test_file, buildFolder+"test.f90")
			#if(self.has_plot): shutil.copyfile(plot_file, buildFolder+"plot.gps")
			#print " done!"

		#copy solver files to build folder
		print "- copying solver(s) to /build..."
		if(self.useDvodeF90):
			shutil.copyfile("solver/dvode_f90_m.f90", buildFolder+"dvode_f90_m.f90")
		else:
			shutil.copyfile("solver/opkdmain.f", buildFolder+"opkdmain.f")
			shutil.copyfile("solver/opkda1.f", buildFolder+"opkda1.f")
			shutil.copyfile("solver/opkda2.f", buildFolder+"opkda2.f")
		#copy non-linear equation solver to build folder
		if(self.useNLEQ): shutil.copyfile("solver/nleq_all.f", buildFolder+"nleq_all.f")

		#copy utility to list the user functions
		fname = "tools/list_user_functions.py"
		if(os.path.exists(fname)):
			shutil.copyfile(fname, buildFolder+"list_user_functions.py")


	#######################################################
	def indent(self):
		buildFolder = self.buildFolder
		print "Indenting...",
		if(self.doIndent):
			if(self.buildCompact):
				indentF90(buildFolder+"krome_all.f90")
			else:
				indentF90(buildFolder+"krome_user_commons.f90")
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
	def ramses_patch2011(self):
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
		self.replacein(pfold+fname,buildFolder+fname,["#KROME_nvar"],\
			["#this must be NDIM+"+str(ramses_offset)+"+"+str(chemCount)], False)

		#move the krome files into the ramses patch folder
		shutil.move(buildFolder+"krome_all.f90", ramsesFolder+"krome_all.f90")
		shutil.move(buildFolder+"krome_user_commons.f90", ramsesFolder+"krome_user_commons.f90")
		shutil.move(buildFolder+"opkda1.f", ramsesFolder+"opkda1.f")
		shutil.move(buildFolder+"opkda2.f", ramsesFolder+"opkda2.f")
		shutil.move(buildFolder+"opkdmain.f", ramsesFolder+"opkdmain.f")


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
                cheminit = "\n"
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
			cheminit += "q(1:nn,ndim+"+ramses_offset+"+"+str(ichem)+")  = "+sdef+"  !"+x.name+"\n"
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
				bkscaleueq += "unoneq("+str(ichem)+") = unoneq("+str(ichem)+")*"+str(x.mass)+"*iscale_d !"+x.name+"\n"
				bkupdateueq += "uold(ind_leaf(i),ndim+"+ramses_offset+"+"+str(ichem)+") = unoneq("+str(ichem)+")\n"
		org = ["#KROME_update_unoneq","#KROME_scale_unoneq","#KROME_backscale_unoneq","#KROME_backupdate_unoneq"]
		new = [updateueq, scaleueq, bkscaleueq, bkupdateueq]
		#replace pragmas (org) with expressions (new)
		self.replacein(pfold+fname, ramsesFolder+fname, org, new)
		indentF90(ramsesFolder+fname)

		#init_flow_fine
		fname = "init_flow_fine.f90"
		init_array = "\n"
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
		
		#read_hydro_params
		fname = "read_hydro_params.f90"
		self.replacein(pfold+fname,ramsesFolder+fname,[],[])
		indentF90(ramsesFolder+fname)

		#Makefile
		fname = "Makefile"
		#note that makefile will be copied in the build folder
		self.replacein(pfold+fname,buildFolder+fname,["#KROME_nvar"],\
			["#this must be NDIM+"+str(ramses_offset)+"+"+str(chemCount)], False)

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
			
			#4. GridKrome.h
			krome_grid_identifya.append("int &"+name+"Num")

		if(speciesCount%3!=0): krome_driver_args += "&"
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
		prags = ["#KROME_args","#KROME_rprec","#KROME_scale","#KROME_minval","#KROME_dom","#KROME_mod"]
		reps = [krome_driver_args,krome_driver_rprec,krome_driver_scale,krome_driver_minval]
		reps += [krome_driver_dom,krome_driver_mod]
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

		#4. replace in GridKrome.h
		fname = "GridKrome.h"
		self.replacein(patchFolder+fname, enzoFolder+fname, ["#KROME_identify"], [krome_grid_identify], False)

		#5. copy others
		fname = "evaluate_tgas.F90"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "InitializeRateData.C"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
		fname = "krome_initab.F90"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
                fname = "notes.txt"
		shutil.copy(patchFolder+fname, enzoFolder+fname)
                fname = "krome_enzo.sh"
		shutil.copy(patchFolder+fname, enzoFolder+fname)

                if(self.useDvodeF90):
                        fname = "kromebuild_dvode.sh"
		        shutil.copy(patchFolder+fname, enzoFolder+"kromebuild.sh")
                else:
		        fname = "kromebuild.sh"
		        shutil.copy(patchFolder+fname, enzoFolder+fname)

		#6. move others
		flist = ["krome_all", "krome_user_commons"]
		for fle in flist:
			shutil.move(buildFolder+fle+".f90", enzoFolder+fle+".F90")
       
                if(self.useDvodeF90):
		        flist = ["dvode_f90_m"]
		        for fle in flist:
                                shutil.move(buildFolder+fle+".f90", enzoFolder+fle+".F90")
                else:
		        flist = ["opkda1", "opkda2", "opkdmain"]
		        for fle in flist:
			        shutil.move(buildFolder+fle+".f", enzoFolder+fle+".F")

		return
	############################################
	#prepare the patches if needed
	def patches(self):
		if(self.doFlash): self.flash_patch()
		if(self.doRamses): self.ramses_patch()
		#if(self.doRamses2011): self.ramses_patch()
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
			print "WARNING: "+str(len(reacts))+" reactions found! Using implicit RHS (option -iRHS)" 
			print "could be more efficient and also allows faster compilation."
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
