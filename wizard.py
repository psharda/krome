#!/usr/bin/python
#!python
from kromelib import *
ys = ["Y","y"]
ns = ["N","n"]
yns = ys + ns
coolList = ["NONE","ATOMIC","H2","HD","DH","DUST","H2GP98","COMPTON","EXPANSION","CIE","CONT","CHEM","DISS","Z"]
coolHelp = ["No cooling","Atomic from Cen 1992","H2 from Glover+2007","HD from Lipovka+2007","Endothermic with thermochemical data","Dust cooling",\
		"H2 from Galli+Palla 1998","Compton","Isothermal expanding gas","Collisional induced","Continuum emission",\
		"Endothermic reactions","Dissociative","Metal"]

heatList = ["NONE","COMPRESS", "PHOTO", "CHEM", "DH", "CR", "PHOTOAV", "PHOTODUST"]
heatHelp = ["No heating","isothermal copmpressional", "photoheating","exothermic reactions", "exothermic from thermochemical data","cosmic rays",\
	"photoheating due to UV photopumping","photoelectric heating from dust"]

print "********************"
print " option wizard!"
print "********************"
rawc = 0
fullcmd = "./krome"


def ask(optNameIn, optDefaultIn, optMessage, optMode="yn", optList=None, helpList=None):
	optName = optNameIn.strip()
	if(optName[0]!="-"): optName = "-"+optName
	iOptMode = optMode.strip().lower()
	optDefault = optDefaultIn.strip()
	optMessage = optMessage[0].upper()+optMessage[1:]
	rwi = ""
	while(rwi==""):
		if(optList!=None):
			for i in range(len(optList)):
				print str(i)+")",optList[i],helpList[i]
		rwi = raw_input("- "+optMessage+" ["+optDefault+"]: ").strip()
		if(rwi==""): rwi = optDefault
		if(iOptMode=="file"):
			if(file_exists(rwi)):
				optOut = optName + "=" + rwi
			else:
				sys.exit("ERROR: file "+rwi+" not found!")
		elif(iOptMode=="yn"):
			optOut = ""
			if(rwi in ys): optOut = optName
		elif(iOptMode=="list"):
			optOut = optName+"="+(",".join([optList[int(x.strip())] for x in rwi.split(",")]))
			if(optOut==(optName+"="+"NONE")): optOut = ""
		elif(iOptMode=="arg"):
			optOut = ""
			if(rwi!="ignore"):
				optOut = optName+" "+rwi
		else:
			sys.exit("ERROR: unknow mode "+iOptMode+" for "+optNameIn+"!")
	if(optName!="-adv"):
		if(optOut==""):
			print "option "+optName+" ignored\n"
		else:
			print "added "+optOut+"\n"
	return optOut+" "

fullcmd = ""

fullcmd += ask("-n", "networks/react_COthin", "path of your chemical network", "arg")
fullcmd += ask("-useN", "y","Use number density (otherwise mass fractions)?", "yn")
fullcmd += ask("-coolFile", "data/coolZ.dat","Do you want to extend cooling with an additional cooling file?", "file")
fullcmd += ask("-cooling", "0","Cooling functions (use numbers above comma separated)?","list",coolList,coolHelp)
fullcmd += ask("-heating", "0","Heating functions (use numbers above comma separated)?","list",heatList,heatHelp)

adv = ask("-adv", "n","do you want to use advanced THERMAL options?","yn")
if(adv.strip()!=""):
	fullcmd += ask("-coolLevels", "ignore", "Maximum number of cooling lines", "arg")
	fullcmd += ask("-coolingQuench", "ignore", "Quench cooling below the indicated Tgas", "arg")
	fullcmd += ask("-H2opacity", "ignore", "Use opacity for H2. Types are OMUKAI or RIPAMONTI", "arg")
	fullcmd += ask("-gamma", "ignore", "Adiabatic index. Allowed: FULL, VIB, ROT, EXACT, REDUCED, or custom expression.", "arg")
	fullcmd += ask("-useThermoToggle", "n", "Include the possibility of switching off the dT/dt calculation.", "yn")
	fullcmd += ask("-useCoolCMBFloor", "n", "Subract dT/dt evaluated at Tcmb from the one evaluated at Tgas. ", "yn")

adv = ask("-adv", "n","do you want to use advanced SOLVER options?","yn")
if(adv.strip()!=""):
	fullcmd += ask("-ATOL", "ignore", "Set asbolute tolerance", "arg")
	fullcmd += ask("-RTOL", "ignore", "Set relative tolerance", "arg")
	fullcmd += ask("-customATOL", "ignore", "Use individual ATOLs from a file", "arg")
	fullcmd += ask("-customRTOL", "ignore", "Use individual RTOLs from a file", "arg")
	fullcmd += ask("-customODE", "ignore", "Use custom ODE from a file", "arg")
	fullcmd += ask("-forceMF21", "n", "Force the option MF21 of the solver", "yn")
	fullcmd += ask("-forceRWORK", "ignore", "Use a custom value for RWORK instead automatic", "arg")
	fullcmd += ask("-maxord", "ignore", "Maximum order of the solver", "arg")
	fullcmd += ask("-useODEConstant", "ignore", "Add a constant to all the ODEs", "arg")
	fullcmd += ask("-useCustomCoe", "ignore", "Use a custom expression for rate coefficients", "arg")
	fullcmd += ask("-useEquilibrium", "n", "When equilibrium is reached stop integration - DANGEROUS!", "yn")
	fullcmd += ask("-useDvodeF90", "n", "Use dvode.f90 instead of DLSODES", "yn")

adv = ask("-adv", "n","do you want to use advanced BUILD options?","yn")
if(adv.strip()!=""):
	fullcmd += ask("-compact", "n", "Use a single module file named krome_all.f90", "yn")
	fullcmd += ask("-clean", "n", "Wipe the build/ folder before new build", "yn")
	fullcmd += ask("-noExample", "n", "Do not replace test.f90 and Makefile", "yn")
	fullcmd += ask("-dry", "n", "Run KROME but do not write files", "yn")
	fullcmd += ask("-sh", "n", "Use short comment header", "yn")
	fullcmd += ask("-project", "ignore", "Set a name for the build folder as build_NAME", "arg")
	fullcmd += ask("-source", "ignore", "Use a different source folder instead of src", "arg")
	fullcmd += ask("-unsafe", "n", "Do not write warning message", "yn")

adv = ask("-adv", "n","do you want to use advanced NETWORK options?","yn")
if(adv.strip()!=""):
	fullcmd += ask("-nuclearMult", "n", "keep into account reactants multeplicity", "yn")
	fullcmd += ask("-skipDup", "n", "skip duplicate reactions", "yn")
	fullcmd += ask("-Tlimit", "ignore", "set upper and lower temperature operators, comma separated. LE,LT,GE,GT.", "arg")
	fullcmd += ask("-useFileIdx", "n", "use the index found in the file (otherwise automatic)", "yn")
	fullcmd += ask("-noTlimits", "n", "do not use temeperature limits", "yn")
	fullcmd += ask("-mergeTlimits", "n", "reactions with same reactants and products use the same index", "yn")


adv = ask("-adv", "n","do you want to use advanced DEBUG options?","yn")
if(adv.strip()!=""):
	fullcmd += ask("-checkConserv", "n", "check conservation during the integration (slow)", "yn")
	fullcmd += ask("-checkReverse", "n", "check if reverse reactions are present in the network", "yn")
	fullcmd += ask("-conserve", "n", "conserve charge and mass", "yn")
	fullcmd += ask("-conserveE", "n", "conserve charge only", "yn")
	fullcmd += ask("-nochargeCheck", "n", "do not check the charge conservation", "yn")
	fullcmd += ask("-noCheck", "n", "do not check the charge and mass conservation", "yn")
	fullcmd += ask("-nomassCheck", "n", "do not check mass conservation", "yn")
	fullcmd += ask("-pedantic", "n", "copy a pedantic Makefile in the build directory", "yn")
	fullcmd += ask("-report", "n", "dump a report after every call of the solver (very slow!)", "yn")
	fullcmd += ask("-ierr", "n", "use add a ierr integer in the call to krome to handle errors if necessary", "yn")
	fullcmd += ask("-checkThermochem", "n", "check if species are in the thermochemistry database", "yn")


print "\nYour call to krome is:"
print "./krome "+fullcmd
sys.exit()




#PATCH
rawc += 1
patch = ""
while(patch==""):
	for i in range(len(patchList)):
		print str(i)+")",patchList[i],patchHelp[i]
	rwi = raw_input(str(rawc)+'. do you want to prepare a patch? [0]: ').strip()
	if(rwi==""): rwi = "0"

if(rwi!="NONE"): fullcmd += " -heating="+heats

print
print fullcmd


