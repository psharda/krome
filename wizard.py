#!/usr/bin/python
#!python
from kromelib import *
import re,sys,time

outFileName = "options.kop"
argv = sys.argv
if(len(argv)>1):
	if("-h" in argv):
		sys.exit("usage: "+argv[0]+" [OUTPUTFILE] [-h]")
	outFileName = argv[1].strip()

ys = ["y","yes"]
ns = ["N","n"]
yns = ys + ns
coolList = ["NONE","ATOMIC","H2","HD","DH","DUST","H2GP98","COMPTON","EXPANSION","CIE","CONT","CHEM","DISS","Z"]
coolHelp = ["No cooling","Atomic from Cen 1992","H2 from Glover+2007","HD from Lipovka+2007","Endothermic with thermochemical data","Dust cooling",\
		"H2 from Galli+Palla 1998","Compton","Isothermal expanding gas","Collisional induced","Continuum emission",\
		"Endothermic reactions","Dissociative","Metal"]

heatList = ["NONE","COMPRESS", "PHOTO", "CHEM", "DH", "CR", "PHOTOAV", "PHOTODUST"]
heatHelp = ["No heating","isothermal compressional", "photoheating","exothermic reactions", "exothermic from thermochemical data","cosmic rays",\
	"photoheating due to UV photopumping","photoelectric heating from dust"]

dustList = ["NONE,""GROWTH","SPUTTER","H2","T"]
dustHelp = ["no processes","growth from gas","sputter from gas collisions","H2 formation on the surface","grain temperature calculation"]

#atomic number dictionary
zdic = ["E", "D", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S",\
	 "CL", "AR", "K", "TI", "CA", "CR", "MN", "FE", "NI"]
zdic = sorted(zdic,key=lambda x:len(x),reverse=True)

print "********************"
print " option wizard!"
print "********************"
rawc = 0
fullcmd = "./krome"

################
_digits = re.compile('\d')
def contains_digits(d):
    return bool(_digits.search(d))

################
def ask(optNameIn, optDefaultIn, optMessage, optMode="yn", optList=None, helpList=None):
	optName = optNameIn.strip()
	if(optName[0]!="-"): optName = "-"+optName
	iOptMode = optMode.strip().lower()
	optDefault = optDefaultIn.strip()
	optMessage = optMessage[0].upper()+optMessage[1:]
	rwi = ""
	#while you don't write something
	while(rwi==""):
		#if list is available print it
		if(optList!=None):
			for i in range(len(optList)):
				print str(i)+")",optList[i],"("+helpList[i]+")"
		#read input from terminal
		rwi = raw_input("- "+optMessage+" ["+optDefault+"]: ").strip()
		#empty input uses default
		if(rwi==""): rwi = optDefault
		#mode file check if file exists
		if(iOptMode=="file"):
			if(file_exists(rwi)):
				optOut = optName + "=" + rwi
			else:
				sys.exit("ERROR: file "+rwi+" not found!")
		#mode yes/no check if input is yes or y
		elif(iOptMode=="yn"):
			optOut = ""
			if(rwi.lower() in ys): optOut = optName
		#mode list extrat data
		elif(iOptMode=="list"):
			#split comma
			optListInt = [int(x.strip()) for x in rwi.split(",")]
			#use list of unique parameters
			uOptListInt = []
			for x in optListInt:
				if(x in uOptListInt): continue
				uOptListInt.append(x)
			optListInt = uOptListInt
			#check if NONE is included in the list and rise an error
			if((len(optListInt)>1) and (0 in optListInt) and (optList[0]=="NONE")):
				sys.exit("ERROR: you can't include NONE and other options!")
			#join the list of the selected options
			optOut = optName+"="+(",".join([optList[x] for x in optListInt]))
			#if NONE no need to include the option
			if(optOut==(optName+"="+"NONE")): optOut = ""
		#mode argument read the argument and append. If "ignore" skip the option
		elif(iOptMode=="arg"):
			optOut = ""
			if(rwi!="ignore"):
				optOut = optName+"="+rwi
		else:
			sys.exit("ERROR: unknow mode "+iOptMode+" for "+optNameIn+"!")
	#if the request is for an advanced option then no verbose output
	if(optName!="-adv"):
		if(optOut==""):
			print "option "+optName+" ignored\n"
		else:
			print "added "+optOut+"\n"
	return optOut

fullcmd = []

fullcmd.append(ask("-n", "networks/react_COthin", "path of your chemical network", "arg"))
fullcmd.append(ask("-useN", "y","Use number density (otherwise mass fractions)?", "yn"))
coolFileOpt = ask("-coolFile", "data/coolZ.dat","Do you want to extend cooling with an additional cooling file?", "file")
fullcmd.append(coolFileOpt)
#add extra cooling to the cooling list from file
if(coolFileOpt.strip()!=""):
	#grab cooling fname
	fname = coolFileOpt.replace("-coolFile=","").strip()
	#open and read
	fh = open(fname,"rb")
	#look for coolants
	for row in fh:
		srow = row.strip()
		if(not("metal:" in srow) and not("coolant:" in srow)): continue
		metal = srow.replace("metal:","").replace("coolant:","").strip()
		mymet = metal.upper()
		is_atom = True
		#if contains digits is a molecule
		mymet.replace("+","").replace("-","")
		if(contains_digits(mymet)):
			is_atom = False
		else:
			repcount = 0 #replacement counts
			#loop on atom dictionary
			for z in zdic:
				if(z in mymet):
					repcount += 1
					mymet = mymet.replace(z,"")
				if(repcount>1): break #save iterations
			#if more than 1 atom replacement is a molecule
			if(repcount>1): is_atom = False
		coolname = metal
		if(is_atom):
			coolname = getRomanName(metal)
		coolList.append(coolname)
		coolHelp.append(coolname+" cooling")

coolOpt = ask("-cooling", "0","Cooling functions (use numbers above comma separated)?","list",coolList,coolHelp)
fullcmd.append(coolOpt)
#if no cooling remove file cooling option (no needed)
if(coolOpt.strip()==""):
	fullcmd = [x for x in fullcmd if x!=coolFileOpt]

fullcmd.append(ask("-heating", "0","Heating functions (use numbers above comma separated)?","list",heatList,heatHelp))

adv = ask("-adv", "n","do you want to use advanced THERMAL options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-coolLevels", "ignore", "Maximum number of cooling lines", "arg"))
	fullcmd.append(ask("-coolingQuench", "ignore", "Quench cooling below the indicated Tgas", "arg"))
	fullcmd.append(ask("-H2opacity", "ignore", "Use opacity for H2. Types are OMUKAI or RIPAMONTI", "arg"))
	fullcmd.append(ask("-gamma", "ignore", "Adiabatic index. Allowed: FULL, VIB, ROT, EXACT, REDUCED, or custom expression.", "arg"))
	fullcmd.append(ask("-useThermoToggle", "n", "Include the possibility of switching off the dT/dt calculation.", "yn"))
	fullcmd.append(ask("-useCoolCMBFloor", "n", "Subract dT/dt evaluated at Tcmb from the one evaluated at Tgas. ", "yn"))

adv = ask("-adv", "n","do you want to use advanced SOLVER options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-ATOL", "ignore", "Set asbolute tolerance", "arg"))
	fullcmd.append(ask("-RTOL", "ignore", "Set relative tolerance", "arg"))
	fullcmd.append(ask("-customATOL", "ignore", "Use individual ATOLs from a file", "arg"))
	fullcmd.append(ask("-customRTOL", "ignore", "Use individual RTOLs from a file", "arg"))
	fullcmd.append(ask("-customODE", "ignore", "Use custom ODE from a file", "arg"))
	fullcmd.append(ask("-forceMF21", "n", "Force the option MF21 of the solver", "yn"))
	fullcmd.append(ask("-forceRWORK", "ignore", "Use a custom value for RWORK instead automatic", "arg"))
	fullcmd.append(ask("-maxord", "ignore", "Maximum order of the solver", "arg"))
	fullcmd.append(ask("-useODEConstant", "ignore", "Add a constant to all the ODEs", "arg"))
	fullcmd.append(ask("-useCustomCoe", "ignore", "Use a custom expression for rate coefficients", "arg"))
	fullcmd.append(ask("-useEquilibrium", "n", "When equilibrium is reached stop integration - DANGEROUS!", "yn"))
	fullcmd.append(ask("-useDvodeF90", "n", "Use dvode.f90 instead of DLSODES", "yn"))

adv = ask("-adv", "n","do you want to use advanced BUILD options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-compact", "n", "Use a single module file named krome_all.f90", "yn"))
	fullcmd.append(ask("-clean", "n", "Wipe the build/ folder before new build", "yn"))
	fullcmd.append(ask("-noExample", "n", "Do not replace test.f90 and Makefile", "yn"))
	fullcmd.append(ask("-dry", "n", "Run KROME but do not write files", "yn"))
	fullcmd.append(ask("-sh", "n", "Use short comment header", "yn"))
	fullcmd.append(ask("-project", "ignore", "Set a name for the build folder as build_NAME", "arg"))
	fullcmd.append(ask("-source", "ignore", "Use a different source folder instead of src", "arg"))
	fullcmd.append(ask("-unsafe", "n", "Do not write warning message", "yn"))

adv = ask("-adv", "n","do you want to use advanced NETWORK options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-nuclearMult", "n", "keep into account reactants multeplicity", "yn"))
	fullcmd.append(ask("-skipDup", "n", "skip duplicate reactions", "yn"))
	fullcmd.append(ask("-Tlimit", "ignore", "set upper and lower temperature operators, comma separated. LE,LT,GE,GT.", "arg"))
	fullcmd.append(ask("-useFileIdx", "n", "use the index found in the file (otherwise automatic)", "yn"))
	fullcmd.append(ask("-noTlimits", "n", "do not use temeperature limits", "yn"))
	fullcmd.append(ask("-mergeTlimits", "n", "reactions with same reactants and products use the same index", "yn"))
	fullcmd.append(ask("-useTabs", "n", "use interpolated tables for the rate coefficents", "yn"))
	fullcmd.append(ask("-reverse", "n", "create reverse reactions using NASA polynomials", "yn"))
	fullcmd.append(ask("-listAutomatics", "n", "list all the automatic species (large output)", "yn"))
	fullcmd.append(ask("-iRHS", "n", "implicit loop for reaction instead of explicit fex", "yn"))
	fullcmd.append(ask("-usePlainIsotopes ", "n", "use kA format instead of [k]A format, e.g. 13CO instead of [13]CO", "yn"))


adv = ask("-adv", "n","do you want to use advanced DEBUG options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-checkConserv", "n", "check conservation during the integration (slow)", "yn"))
	fullcmd.append(ask("-checkReverse", "n", "check if reverse reactions are present in the network", "yn"))
	fullcmd.append(ask("-conserve", "n", "conserve charge and mass", "yn"))
	fullcmd.append(ask("-conserveE", "n", "conserve charge only", "yn"))
	fullcmd.append(ask("-nochargeCheck", "n", "do not check the charge conservation", "yn"))
	fullcmd.append(ask("-noCheck", "n", "do not check the charge and mass conservation", "yn"))
	fullcmd.append(ask("-nomassCheck", "n", "do not check mass conservation", "yn"))
	fullcmd.append(ask("-pedantic", "n", "copy a pedantic Makefile in the build directory", "yn"))
	fullcmd.append(ask("-report", "n", "dump a report after every call of the solver (very slow!)", "yn"))
	fullcmd.append(ask("-ierr", "n", "use add a ierr integer in the call to krome to handle errors if necessary", "yn"))
	fullcmd.append(ask("-checkThermochem", "n", "check if species are in the thermochemistry database", "yn"))

adv = ask("-adv", "n","do you want to use DUST options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-dust", "ignore", "set the number of dust bins and the type(s), e.g. 10,C,Si", "arg"))
	fullcmd.append(ask("-dustOptions", 0, "include dust-related processes", "list",dustList,dustHelp))

adv = ask("-adv", "n","do you want to use PHOTOCHEMISTRY options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-photoBins", "ignore", "number of photon energy bins employed", "arg"))
	fullcmd.append(ask("-shielding", "ignore", "type of shielding employed for photochemistry, WG11 or DB96","arg"))
	fullcmd.append(ask("-usePhotoInduced", "n", "use photon-induced transitions in the cooling","yn"))
	fullcmd.append(ask("-usePhotoOpacity", "n", "compute photo reactions at runtime using chemical-dependent opacity","yn"))

fullcmd.append(ask("-stars", "n", "enable chemistry for star thermochemistry","yn"))

adv = ask("-adv", "n","do you want to use PATCH options?","yn")
if(adv.strip()!=""):
	fullcmd.append(ask("-C", "n", "enable patch for generic C codes (no needed for ENZO)","yn"))
	fullcmd.append(ask("-enzo", "n", "build a patch for ENZO","yn"))
	fullcmd.append(ask("-flash", "n", "build a patch for FLASH","yn"))
	fullcmd.append(ask("-ramses", "n", "build a patch for RAMSES","yn"))
	fullcmd.append(ask("-ramsesTH", "n", "build a patch for private version of RAMSES","yn"))
	fullcmd.append(ask("-ramsesOffset", "ignore", "offset for RAMSES NVAR","arg"))
	

fullcmd = [x for x in fullcmd if(x!="")]


print "\nYour call to krome is:"
print " ./krome "+(" ".join(fullcmd))

print "or use the option file"
print " "+outFileName
fh = open(outFileName,"w")
fh.write("#This option file has been automatically generated\n")
fh.write("#"+time.asctime(time.localtime(time.time()))+"\n\n")
for opt in fullcmd:
	fh.write(opt+"\n\n")
	


