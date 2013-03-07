import sys,shutil,os,glob
from kromelib import *

#check number the of arguments
if(len(sys.argv)<2):
	die(get_usage())


#you can select only one -forceMF
if(("-forceMF=222" in sys.argv) and ("-forceMF=21" in sys.argv)):
	die("ERROR: options -forceMF=222 and -forceMF=21 are mutually exclusive: choose one.")
print

#set defaults
solver_MF = 222
force_rwork = useHeating = doReport = checkConserv = useFileIdx = buildCompact = False
use_implicit_RHS = use_photons = useTabs = useDvodeF90 = useTopology = useFlux = False
useCoolingCEN = useCoolingH2 = useCoolingH2GP98 = useCoolingHD = use_cooling = False
createReverse = useCustomCoe = useODEConstant = cleanBuild = usePlainIsotopes = useDust = False
useX = has_plot = True
test_name = "default"
is_test = False
TlimitOpLow = "GE"
TlimitOpHigh = "LT"
customCoeFunction = "[CUSTOM COE FUNCTION NOT SET!]"
buildFolder = "build/"
TminAuto = 1e99
TmaxAuto = -1e99
dustArraySize = dustTypesSize = 0
#select test name
for arg in sys.argv:
	if("test=" in arg):
		is_test = True
		test_name = (arg.strip().replace("-test=",""))
		print "Reading option -test (test="+test_name+")"
		if(test_name=="planet"):
			[sys.argv.append(x) for x in ["-useN"]]
			filename = "networks/react_planet"
		elif(test_name=="WH2008"):
			[sys.argv.append(x) for x in ["-useN","-iRHS"]]
			filename = "networks/react_WH2008"
		elif(test_name=="fake"):
			[sys.argv.append(x) for x in ["-forceRWORK=200"]]
			filename = "networks/react_fake"
		elif(test_name=="topology"):
			[sys.argv.append(x) for x in ["-useN","-iRHS","-useTopology"]]
			filename = "networks/react_WH2008"
		elif(test_name=="flux"):
			[sys.argv.append(x) for x in ["-useN","-iRHS","-useFlux"]]
			filename = "networks/react_WH2008"
		elif(test_name=="cooling"):
			[sys.argv.append(x) for x in ["-useN","-cooling=ALL"]]
			filename = "networks/react_enzo"
		elif(test_name=="shock1Dcool"):
			test_name = "shock1D"
			[sys.argv.append(x) for x in ["-cooling=CEN,H2,HD"]]
			filename = "networks/react_enzo"
		elif(test_name=="shock1D"):
			filename = "networks/react_enzo"
		elif(test_name=="dust"):
			has_plot= False
			[sys.argv.append(x) for x in ["-dust=10,C,Si","-useN"]]
			filename = "networks/react_enzo"
		else:
			print "ERROR: test \""+test_name+"\" not present!"
			sys.exit()
		break
#get filename
if(not(is_test)): filename = sys.argv[1]

#set implicit RHS
if("-iRHS" in sys.argv):
	use_implicit_RHS = True
	solver_MF = 222
	print "Reading option -iRHS"
#force MF=21
if("-forceMF=21" in sys.argv):
	solver_MF = 21
	print "Reading option -forceMF=21"
#force MF=222
if("-forceMF=222" in sys.argv):
	solver_MF = 222
	print "Reading option -forceMF=222"
#use numeric density instead of fractions as input
if("-useN" in sys.argv):
	useX = False
	print "Reading option -useN"
#force to use photons
if("-usePhot" in sys.argv):
	use_photons = True
	print "Reading option -usePhot"
#use rate tables
if("-useTabs" in sys.argv):
	useTabs = True
	print "Reading option -useTabs"
#use f90 solver
if("-useDvodeF90" in sys.argv):
	useDvodeF90 = True
	print "Reading option -useDvodeF90"
#use topology reduction
if("-useTopology" in sys.argv):
	useTopology = True
	print "Reading option -useTopology"
#use topology reduction
if("-useFlux" in sys.argv):
	useFlux = True
	print "Reading option -useFlux"
#use heating
if("-useHeating" in sys.argv):
	useHeating = True
	print "Reading option -useHeating"
#do report
if("-report" in sys.argv):
	doReport = True
	print "Reading option -report"
#check mass conservation
if("-checkConserv" in sys.argv):
	checkConserv = True
	print "Reading option -checkConserv"
#use reaction indexes in reaction file
if("-useFileIdx" in sys.argv):
	useFileIdx = True
	print "Reading option -useFileIdx"
#write a single compact file krome_all.f90
if("-compact" in sys.argv):
	buildCompact = True
	print "Reading option -compact"
#perform a clean build
if("-clean" in sys.argv):
	cleanBuild = True
	print "Reading option -clean"
#build isotopes automatically
if("-usePlainIsotopes" in sys.argv):
	usePlainIsotopes = True
	print "Reading option -usePlainIsotopes"
#determine reverse function to reverse reactions
for arg in sys.argv:
	if("reverse=" in arg):
		reverse_function = (arg.strip().replace("-reverse=",""))
		createReverse = True
		print "Reading option -reverse (function="+reverse_function+")"
		break
#determine Tgas limit operators
for arg in sys.argv:
	if("Tlimit=" in arg):
		myTlimit = (arg.strip().replace("-Tlimit=",""))
		myTlimit = myTlimit.replace("[","").replace("]","").split(",")
		TlimitOpHigh = myTlimit[1].strip().upper()
		TlimitOpLow = myTlimit[0].strip().upper()
		allOps = ["LE","LT","GE","GT"]
		if(not(TlimitOpLow in allOps) or not(TlimitOpHigh in allOps)):
			die("ERROR: on -Tlimit operators must be one of the followings: "+(", ".join(allOps)))
		print "Reading option -Tlimit (Low="+TlimitOpLow+", High="+TlimitOpHigh+")"
		break
#determine cooling types
for arg in sys.argv:
	if("cooling=" in arg):
		myCools = arg.strip().replace("-cooling=","").upper().split(",")
		allCools = ["CEN","H2","H2GP98","HD","ALL"]
		for coo in myCools: 
			if(not(coo in allCools)): 
				die("ERROR: Cooling \""+coo+"\" is unknown!\nAvailable coolings are: "+(", ".join(allCools)))

		if("ALL" in myCools): myCools = allCools
		if("CEN" in myCools): useCoolingCEN = True
		if("H2" in myCools): useCoolingH2 = True
		if("H2GP98" in myCools): useCoolingH2GP98 = True
		if("HD" in myCools): useCoolingHD = True
		use_cooling = True

		print "Reading option -cooling ("+(",".join(myCools))+")"
		break

#force rwork size
for arg in sys.argv:
	if("testFile=" in arg):
		myrwork = (arg.strip().replace("-testFile=",""))
		is_test = True
		test_name = myrwork
		print "Reading option -testFile (file="+str(myrwork)+")"
		break

#force rwork size
for arg in sys.argv:
	if("forceRWORK=" in arg):
		myrwork = (arg.strip().replace("-forceRWORK=",""))
		force_rwork = True
		print "Reading option -forceRWORK (RWORK="+str(myrwork)+")"
		break
#use custom function for coefficient instead of coe_tab()
for arg in sys.argv:
	if("useCustomCoe=" in arg):
		customCoeFunction = (arg.strip().replace("-useCustomCoe=",""))
		useCustomCoe = True
		print "Reading option -useCustomCoe (Expression="+str(customCoeFunction)+")"
		break
#use function to append after each ODE
for arg in sys.argv:
	if("useODEConstant=" in arg):
		ODEConstant = (arg.strip().replace("-useODEConstant=",""))
		useODEConstant = True
		print "Reading option -useODEConstant (Constant="+str(ODEConstant)+")"
		break
#use function to append after each ODE
for arg in sys.argv:
	if("dust=" in arg):
		dustopt = (arg.strip().replace("-dust=",""))
		adust = dustopt.split(",")
		useDust = True
		if(len(adust)<2): die("ERROR: you must specify dust size and type(s), e.g. -dust=20,C,Si")
		dustArraySize = int(adust[0])
		dustTypes = adust[1:]
		dustTypesSize = len(dustTypes)
		print "Reading option -dust (size="+str(dustArraySize)+", type(s)="+(",".join(dustTypes))+")"
		break
#use function to append after each ODE
for arg in sys.argv:
	if("project=" in arg):
		projectName = (arg.strip().replace("-project=",""))
		print "Reading option -project (name="+str(projectName)+")"
		buildFolder = "build_"+projectName+"/"
		fout = open(projectName+".kpj","w")
		fout.write((" ".join(sys.argv)))
		fout.close()
		break

#show help
if("-help" in sys.argv or "-h" in sys.argv or "--help" in sys.argv):
	die(get_usage())
print

use_RHS_variable = False

separator = "," #separator character
krome = "" #krome global variable prefix

#chech if reactions file exists
try:
	with open(filename) as f: pass
except IOError as e:
	print "ERROR: Reaction file \""+filename+"\" doesn't exist!"
	sys.exit()


solver_moss = int(solver_MF/100)
solver_meth = int(solver_MF/10) % 10
solver_miter = solver_MF % 10

print "solver info:"
print " MF:",solver_MF
print " MOSS+METH+MITER:","+".join([str(x) for x in [solver_moss,solver_meth,solver_miter]])
print

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
if(len(atoms_iso)!=len(atoms_p)): die("ERROR: in building isotopes the length of the atoms array and the number of protons array mismatch!")
for aiso in atoms_iso:
	protons = atoms_p[atoms_iso.index(aiso)] #get proton numbers
	for i in range(protons,80):
		iso_name = "["+str(i)+aiso+"]"
		if(usePlainIsotopes): iso_name = str(i)+aiso
		mass_dic[iso_name] = protons*(me+mp) + ((i-protons)*mn)

##################################
##################################
##################################

mass_dic = dict([[k.upper(),v] for (k,v) in mass_dic.iteritems()])
atoms = sorted(mass_dic, key = lambda x: len(x),reverse=True)

print "Reading from file \""+filename+"\"..."
spec_names = [] #string
specs = []
reacts = []
idx_list = [] #store reaction index in case of -useFileIdx
fh = open(filename,"rb")
rcount = 0 #count reactions
found_one = False #flag to control if at least one reaction has been found
unmatch_idx = False #controls if the reaction index in the file match the sequential index
for row in fh:
	srow = row.strip()
	arow = srow.split(separator,10) #split only 10 elements
	if(srow[0]=="#"): continue #looks for comment line
	if(len(arow)!=11): continue #check line format
	found_one = True
	rcount += 1
	myrea = reaction() #create objec reaction

	if(useFileIdx):
		reaction_idx = int(arow[0]) #index of the reaction (from file)
		if(reaction_idx<1): die("ERROR: reaction index must be > 0!\n Check your reaction file!")
		if(reaction_idx in idx_list): 
			die("ERROR: reaction index "+str(reaction_idx)+" already present!\n Check your reaction file!")
		myrea.idx = reaction_idx
		rcount = reaction_idx
	else:
		myrea.idx = rcount
		if(rcount!=int(arow[0])): unmatch_idx = True
	reactants = arow[1:4] #get reactants
	products = arow[4:8] #get products
			
	myrea.Tmin = format_double(arow[8]) #get Tmin
	myrea.Tmax = format_double(arow[9]) #get Tmax
	TminAuto = min(float(arow[8].lower().replace("d","e")), TminAuto)
	TmaxAuto = max(float(arow[9].lower().replace("d","e")), TmaxAuto)
	myrea.krate = arow[10] #get reaction rate written in F90 style
	if(useCustomCoe): myrea.krate = "0.d0" #when custom function is used standard coefficient are set to zero
	#loop over reactants to grep molecules
	for r in reactants:
		if(r.strip()=="g" and not(use_photons)): continue
		if(r.strip()!=""):
			mol = parser(r,mass_dic,atoms)
			if(not(mol.name in spec_names)):
				spec_names.append(mol.name)
				specs.append(mol)
			mol.idx = spec_names.index(mol.name) + 1
			myrea.reactants.append(mol) #add molecule object to reactants
	#loop over prodcuts to grep molecules
	for p in products:
		if(p.strip()=="g" and not(use_photons)): continue
		if(p.strip()!=""):
			mol = parser(p,mass_dic,atoms)
			if(not(mol.name in spec_names)):
				spec_names.append(mol.name)
				specs.append(mol)
			mol.idx = spec_names.index(mol.name) + 1
			myrea.products.append(mol) #add molecule object to products

	myrea.build_verbatim() #build reaction as string (e.g. A+B->C)
	myrea.reactants = sorted(myrea.reactants, key=lambda r:r.idx) #sort reactants
	myrea.products = sorted(myrea.products, key=lambda p:p.idx) #sort products
	myrea.build_RHS() #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
	myrea.check() #check mass and charge conservation
	reacts.append(myrea)

#check file format
if(not(found_one)):
	die("ERROR: no reactions found in file \""+filename+"\"")
if(unmatch_idx):
	print "WARNING: index in \""+filename+"\" are not sequential!"

print "done!"


#do reverse reaction if needed
if(createReverse):
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
			myrev.build_verbatim() #build reaction as string (e.g. A+B->C)
			myrev.verbatim += " (REV)" #append REV label to reaction string
			myreactants = myrev.reactants #get list of reactants
			myproducts = myrev.products #get list of products
			myrev_function = reverse_function #make a copy of the function
			#loop over reactants and replace template name nRi with local name [e.g. nR1->n(idx_H)]
			for i in range(1,10):
				ops = ["*","/","+","-"]
				#replace reactants token (remove dummy)
				if(i<len(myreactants)+1):
					myr = myreactants[i-1]
					myrev_function = myrev_function.replace("nR"+str(i),"n("+myr.fidx+")")
				else:
					#loop on operators to remove token properly
					for op in ops:
						myrev_function = myrev_function.replace(op+"nR"+str(i),"")
				
				#replace products token (remove dummy)
				if(i<len(myproducts)+1):
					myp = myproducts[i-1]
					myrev_function = myrev_function.replace("nP"+str(i),"n("+myp.fidx+")")
				else:
					#loop on operators to remove token properly
					for op in ops:
						myrev_function = myrev_function.replace(op+"nP"+str(i),"")

			myrev.krate = myrea.krate + myrev_function #get reaction rate written in F90 style
			myrev.build_RHS() #build RHS in F90 format (e.g. k(2)*n(10)*n(8) )
			myrev.check() #check mass and charge conservation
			reacts.append(myrev)
	print "Inverse reaction added: "+str(count_reverse)


#add dust to problem
if(useDust):
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
				mymol = parser(dType,mass_dic,atoms)
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

	print "Dust added:",dustArraySize*dustTypesSize

#look for photons and CR
has_g = has_CR = False

for mol in specs:
	if(mol.name=="g"):
		has_g = True
	if(mol.name=="CR"):
		has_CR = True

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
	photon_species = mymol #store object
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
	photon_species = mymol #store object
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
Tgas_species = mymol #store object

#append dummy as species named dummy
mymol = molec()
mymol.name = "dummy"
mymol.charge = 0
mymol.mass = 0.e0
mymol.ename = "dummy"
mymol.fidx = "idx_dummy"
mymol.idx = len(specs)+1
specs.append(mymol)
dummy = mymol #store object

#print found values
print
print "ODEs needed:", len(specs)
print "Reactions found:", len(reacts)
print "Species found:", len(specs)-4-dustArraySize*dustTypesSize

#dump species to log file
fout = open("species.log","w")
idx = 0
for mol in specs:
	idx += 1
	fout.write(str(idx)+"\t"+mol.name+"\t"+mol.fidx+"\n")
fout.close()

print
#include dust in ODE array partition
sdust = ""
if(useDust):
	#include dust by types
	for dType in dustTypes:
		sdust += " + ["+str(dustArraySize)+" "+dType+"-dust]"
sODEpart = "ODE partition: [" + str(len(specs)-4-dustArraySize*dustTypesSize) + " atom/mols] "+sdust+" + [1 CR] + [1 PHOT] + [1 Tgas] + [1 dummy] = "+str(len(specs))+" ODEs"
print sODEpart
#if species are few print list of species
if(len(specs)<40): print "ODEs list: "+(", ".join([x.name for x in specs]))
print

print "Temperature range found min/max (K):", TminAuto, "/",TmaxAuto

#create explicit differnetials
dns = ["dn("+str(sp.idx)+") = 0.d0" for sp in specs] #initialize
for rea in reacts:
	for r in rea.reactants:
		dns[r.idx-1] = dns[r.idx-1].replace(" = 0.d0"," =")
		dns[r.idx-1] += " -"+rea.RHS
	for p in rea.products:
		dns[p.idx-1] = dns[p.idx-1].replace(" = 0.d0"," =")
		dns[p.idx-1] += " +"+rea.RHS

#create implicit RHS arrays
arr_rr = [[] for i in range(3)]
arr_pp = [[] for i in range(4)]
for rea in reacts:
	for i in range(3):
		if(i<len(rea.reactants)): arr_rr[i].append(rea.reactants[i].idx)
		else: arr_rr[i].append(dummy.idx)
	for i in range(4):
		if(i<len(rea.products)): arr_pp[i].append(rea.products[i].idx)
		else: arr_pp[i].append(dummy.idx)
implicit_arrays = ""
for i in range(len(arr_rr)):
	implicit_arrays += "arr_r"+str(i+1)+"(:) = (/"+(",".join([str(x) for x in arr_rr[i]]))+"/)\n"
for i in range(len(arr_pp)):
	implicit_arrays += "arr_p"+str(i+1)+"(:) = (/"+(",".join([str(x) for x in arr_pp[i]]))+"/)\n"

#wrap RHS (e.g. knn+knn=2*knn)
dnw = []
for dn in dns:
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
		if((i+1) % 4 == 0): dns += "&\n" #break long lines
		if("-" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("-"," -"+str(RHSc[i])+".d0*")
		if("+" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("+"," +"+str(RHSc[i])+".d0*")
		if(RHSc[i]==1): dns+= " "+RHSs[i]
	dnw.append(dns)

#create explicit JAC
neq = len(specs)
jac = [["pdj("+str(j+1)+") = 0.d0" for i in range(neq)] for j in range(neq)]
jsparse = [[0 for i in range(neq)] for j in range(neq)]
for rea in reacts:
	for ri in range(len(rea.reactants)):
		r1 = rea.reactants[ri]
		sjac = "k("+str(rea.idx)+")"
		for rj in range(len(rea.reactants)):
			r2 = rea.reactants[rj]
			if(ri!=rj):
				sjac += "*n("+str(r2.idx)+")"
		for rr in rea.reactants:
			jac[rr.idx-1][r1.idx-1] = jac[rr.idx-1][r1.idx-1].replace(" = 0.d0"," =")
			jac[rr.idx-1][r1.idx-1] += " -"+sjac
			jsparse[rr.idx-1][r1.idx-1] = 1
		for pp in rea.products:
			jac[pp.idx-1][r1.idx-1] = jac[pp.idx-1][r1.idx-1].replace(" = 0.d0"," =")
			jac[pp.idx-1][r1.idx-1] += " +"+sjac
			jsparse[pp.idx-1][r1.idx-1] = 1

for i in range(neq):
	jsparse[i][Tgas_species.idx-1] = 1

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
print "Jacobian non-zero elements:",nnz,"over",neq*neq,"("+str(pnnz)+"% of total elements, sparsity = "+str(100.-pnnz)+"%)"

#if MF=222 no need for sparsity structure arrays
if(solver_MF == 222):
	iaf = jaf = ""
	
#estimate size of RWORK array (see DLSODES manual)
neq = len(specs) #number of equations
nrea = len(reacts) #number of reactions
lenrat = 2
nnz_estimate = pow(neq,2) - nrea #estimate non-zero elements (inaccurate)
nnz = max(nnz, nnz_estimate) #estimate can be more realistic than nnz calculation 
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

if(force_rwork):
	lrw = myrwork
#lrw = int(20+(2+0.5)*nnz + (11+4.5)*neq) #RWORK size

print "LWM:",lwm,"LRW:",lrw

#create build folder if not exists
if(not(os.path.exists(buildFolder))): 
	os.mkdir(buildFolder)
	print "Created "+buildFolder

if(cleanBuild):
	clear_dir(buildFolder) #clear the build directory
	print "Deleted all contents in "+buildFolder
else:
	underFiles = ["commons","constants","cooling","dust","ode","reduction","subs","tabs","user",]
	delFiles = [buildFolder+"krome_"+x+".f90" for x in underFiles] + [buildFolder+x for x in ["krome.f90","opkda1.f","opkda2.f"]]
	delFiles += glob.glob(buildFolder+"*~") + glob.glob(buildFolder+"*.mod") 
	delFiles += glob.glob(buildFolder+"*_genmod.f90") + glob.glob(buildFolder+"*.i90")
	for dfile in delFiles:
		print "deleting "+dfile



print
print "Prepearing files in /build..."

#build all the modules in a single file krome_all.f90
if(buildCompact):
	fout = open(buildFolder+"krome_all.f90","w")

#*********COMMONS****************
#write parameters in krome_commons.f90
print "- writing krome_commons.f90...",
fh = open("src/krome_commons.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_commons.f90","w")
skip = False
for row in fh:

	if(row.strip() == "#KROME_species_index"):
		for x in specs:
			fout.write("\tinteger,parameter::" + x.fidx + "=" + str(x.idx) + "\n")
	elif(row.strip() == "#KROME_parameters"):
			fout.write("\tinteger,parameter::nrea=" + str(len(reacts)) + "\n")
			fout.write("\tinteger,parameter::nmols=" + str(len(specs)-4-dustArraySize*dustTypesSize) + "\n")
			fout.write("\tinteger,parameter::nspec=" + str(len(specs)) + "\n")
			fout.write("\tinteger,parameter::ndust=" + str(dustArraySize*dustTypesSize) + "\n")
			fout.write("\tinteger,parameter::ndustTypes=" + str(dustTypesSize) + "\n")
	elif(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	else:
		if(row[0]!="#"): fout.write(row)
if(not(buildCompact)):
	fout.close()
print "done!"

#********* CONSTANTS ****************
fh = open("src/krome_constants.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_constants.f90","w")
for row in fh:
	if(row[0]!="#"): fout.write(row)
if(not(buildCompact)):
	fout.close()

#*********USER COMMONS****************
#write parameters in krome_user_commons.f90
if(not(file_exists(buildFolder+"krome_user_commons.f90"))):
	print "- writing krome_user_commons.f90...",

	fh = open("src/krome_user_commons.f90")
	if(not(buildCompact)):
		fout = open(buildFolder+"krome_user_commons.f90","w")

	for row in fh:
		if(row.strip() == "#KROME_header"):
			fout.write(get_licence_header())
		else:
			if(row[0]!="#"): fout.write(row)
	if(not(buildCompact)):
		fout.close()
	print "done!"


#*********SUBS****************
#write parameters in krome_subs.f90
print "- writing krome_subs.f90...",
fh = open("src/krome_subs.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_subs.f90","w")
for row in fh:
	if(row.strip() == "#KROME_krates"):
		for x in reacts:
			sTlimit = "if(Tgas."+TlimitOpLow+"."+x.Tmin+" .and. Tgas."+TlimitOpHigh+"."+x.Tmax+")"
			kstr = "\t" + sTlimit + " k("+str(x.idx)+") = " + x.krate + " !" + x.verbatim
			kstr = truncF90(kstr, 60,"*")
			fout.write(truncF90(kstr, 60,"/")+"\n")
	elif(row.strip() == "#KROME_implicit_arrays"):
		fout.write(truncF90(implicit_arrays,60,","))
	elif(row.strip() == "#KROME_masses"):
		for x in specs:
			massrow = "\tget_mass("+str(x.idx)+") = " + str(x.mass).replace("e","d") + "\t!" + x.name + "\n"
			fout.write(massrow.replace("0.0","0.d0"))
	elif(row.strip() == "#KROME_names"):
		for x in specs:
			fout.write("\tget_names("+str(x.idx)+") = \"" + x.name + "\"\n")
	elif(row.strip() == "#KROME_charges"):
		for x in specs:
			fout.write("\tget_charges("+str(x.idx)+") = " + str(x.charge) + ".d0 \t!" + x.name + "\n")
	elif(row.strip() == "#KROME_reaction_names"):
		for x in reacts:
			kstr = "\tget_rnames("+str(x.idx)+") = \"" + x.verbatim +"\""
			fout.write(kstr+"\n")

	elif(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	else:
		fout.write(row)
if(not(buildCompact)):
	fout.close()
print "done!"


#********* TABS ****************
print "- writing krome_tabs.f90...",
fh = open("src/krome_tabs.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_tabs.f90","w")

skip = False
for row in fh:
	if(row.strip() == "#IFKROME_useCustomCoe" and not(useCustomCoe)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_useTabs" and not(useTabs)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_useStandardCoe" and ((useCustomCoe) or (useTabs))): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(skip): continue
	row = row.replace("#KROME_logTlow", "ktab_logTlow = log10(max("+str(TminAuto)+",2.73d0))")
	row = row.replace("#KROME_logTup", "ktab_logTup = log10("+str(TmaxAuto)+")")
	if(useCustomCoe): row = row.replace("#KROMEREPLACE_customCoeFunction", customCoeFunction)

	if(row[0]!="#"): fout.write(row)

if(not(buildCompact)):
	fout.close()
print "done!"
#********* DUST ****************
print "- writing krome_dust.f90...",
fh = open("src/krome_dust.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_dust.f90","w")

skip = False
for row in fh:
	if(row.strip() == "#IFKROME_useDust" and not(useDust)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(skip): continue
	
	dustPartnerIdx = ""
	
	if(useDust):
		itype = 0 #dust type index
		#loop in dust types
		for dType in dustTypes:
			itype += 1 #increase index
			dustPartnerIdx += "krome_dust_partner_idx("+str(itype)+") = idx_"+dType+"\n"
	row = row.replace("#KROME_dustParnerIndex", dustPartnerIdx)

	if(row[0]!="#"): fout.write(row)

if(not(buildCompact)):
	fout.close()
print "done!"

#*********COOLING****************
#write header in krome_cooling.f90
print "- writing krome_cooling.f90...",
fh = open("src/krome_cooling.f90")

if(not(buildCompact)):
	fout = open(buildFolder+"krome_cooling.f90","w")
skip = False
for row in fh:
	if(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	else:
		if(row.strip() == "#IFKROME_useHeating" and not(useHeating)): skip = True
		if(row.strip() == "#ENDIFKROME"): skip = False

		if(row.strip() == "#IFKROME_useCoolingCEN" and not(useCoolingCEN)): skip = True
		if(row.strip() == "#ENDIFKROME"): skip = False

		if(row.strip() == "#IFKROME_useCoolingH2" and not(useCoolingH2)): skip = True
		if(row.strip() == "#ENDIFKROME"): skip = False

		if(row.strip() == "#IFKROME_useCoolingH2GP" and not(useCoolingH2GP98)): skip = True
		if(row.strip() == "#ENDIFKROME"): skip = False

		if(row.strip() == "#IFKROME_useCoolingHD" and not(useCoolingHD)): skip = True
		if(row.strip() == "#ENDIFKROME"): skip = False

		if(skip): continue
		if(row[0]!="#"): fout.write(row)

if(not(buildCompact)):
	fout.close()
print "done!"


#*********ODE****************
#write parameters in krome_ode.f90
print "- writing krome_ode.f90...",

fh = open("src/krome_ode.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_ode.f90","w")
skip = False
for row in fh:
	if(row.strip() == "#IFKROME_use_cooling" and not(use_cooling)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_report" and not(doReport)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(skip): continue
	
	if(row.strip() == "#KROME_ODE"):
		if(use_implicit_RHS):
			fout.write(get_implicit_ode(useTopology or useFlux)+"\n")
		else:
			for x in dnw:
				fout.write("\t" + x + "\n")
	elif(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	elif(row.strip() == "#KROME_implicit_variables"):
		fout.write("real*8::rr\n")
		fout.write("integer::i,r1,r2,r3,p1,p2,p3,p4\n")
	elif(row.strip() == "#KROME_odeConstant" and useODEConstant):
		fout.write("dn(:) = dn(:) "+ODEConstant+" \n")
	elif(row.strip() == "#KROME_odeDust" and useODEConstant):
		fout.write("dn(:) =  krome_dust \n")
	elif(row.strip() == "#KROME_JAC_PD"):
		if(solver_MF==222):
			fout.write("\n")
		else:
			#J(i,j) = df_i/dx_j
			for i in range(neq):
				if(i==0): fout.write("if(j=="+str(i+1)+") then\n")
				if(i>0): fout.write("elseif(j=="+str(i+1)+") then\n")
				if(i!=Tgas_species.idx-1):
					spdj = ""
					has_pdj = False
					for j in range(neq):
						if(jsparse[j][i]==1):
							has_pdj = True
							spdj += ("\t" + jac[j][i] + "\n")
					if(has_pdj): fout.write("k(:) = coe_tab(Tgas)\n")
					fout.write(spdj)
				else:
					jacT = """!rough estimate for dT/dx_i
						do i=1,neq-2
							if(n(i)-nold(i).ne.0.d0) then
								pdj(i) = (n(idx_Tgas) - nold(idx_Tgas))/(n(i) - nold(i))
							end if
						end do"""
					fout.write("\t" + jacT.replace("\t","") + "\n")
			fout.write("end if\n")
	else:
		srow = row.strip()
		if(len(srow)>0):
			if(srow[0]!="#"): fout.write(row)
		else:
			fout.write(row)
if(not(buildCompact)):
	fout.close()
print "done!"

#*********USER****************
#write parameters in krome_user.f90
print "- writing krome_user.f90...",

fh = open("src/krome_user.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_user.f90","w")
for row in fh:

	if(row.strip() == "#IFKROME_use_cooling" and not(use_cooling)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False
	if(skip): continue

	if(row.strip() == "#KROME_species"):
		for x in specs:
			fout.write("\tinteger,parameter::" + "KROME_"+x.fidx + " = " + str(x.idx) +"\t!"+x.name+"\n")
	elif(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	else:
		srow = row.strip()
		if(len(srow)>0):
			if(srow[0]!="#"): fout.write(row)
		else:
			fout.write(row)
if(not(buildCompact)):
	fout.close()
print "done!"

#********* REDUCTION ****************
fh = open("src/krome_reduction.f90")
if(not(buildCompact)):
	fout = open(buildFolder+"krome_reduction.f90","w")
for row in fh:
	if(row[0]!="#"): fout.write(row)
if(not(buildCompact)):
	fout.close()

print "done!"


#*********MAIN****************
#write WORKS arrays and IAC/JAC in krome.f90
print "- writing krome.f90...",
if(useDvodeF90):
	fh = open("src/kromeF90.f90")
else:
	fh = open("src/krome.f90")

if(not(buildCompact)):
	fout = open(buildFolder+"krome.f90","w")
skip = False
for row in fh:
	if(row.strip() == "#IFKROME_useX" and not(useX)): skip = True
	if(row.strip() == "#ELSEKROME" and not(useX)): skip = False
	if(row.strip() == "#ELSEKROME" and useX): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_useTabs" and not(useTabs)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_useFlux" and not(useFlux)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_report" and not(doReport)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_useTopology" and not(useTopology)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(row.strip() == "#IFKROME_check_mass_conservation" and not(checkConserv)): skip = True
	if(row.strip() == "#ENDIFKROME"): skip = False

	if(useDust):
		row = row.replace("#KROME_dust_arguments",",xdust")
	else:
		row = row.replace("#KROME_dust_arguments","")

	if(skip): continue

	if(row.strip() == "#KROME_header"):
		fout.write(get_licence_header())
	elif(row.strip() == "#KROME_rwork_array"):
		fout.write("\treal*8::rwork("+str(lrw)+")\n")
	elif(row.strip() == "#KROME_iwork_array"):
		fout.write("\tinteger::iwork("+str(30+len(ia)+len(ja))+")\n")
	elif(row.strip() == "#KROME_init_IAC"):
		fout.write("\t"+iaf+"\n")
	elif(row.strip() == "#KROME_init_JAC"):
		fout.write("\t"+jaf+"\n")
	elif(row.strip() == "#KROME_MF"):
		if(solver_MF==21):
			fout.write("\tMF=21\n")
		elif(solver_MF==222):
			fout.write("\tMF=222\n")
	else:
		if(row[0]!="#"): fout.write(row)
if(not(buildCompact)):
	fout.close()
print "done!"


if(buildCompact): fout.close()

#*********REPORT.gps****************
#write gnuplot script to plot abundances evolution
if(doReport):
	fout = open(buildFolder+"report.gps","w")
	fout.write("#plot KROME report\n")
	fout.write("set logscale\n")
	fout.write("set xlabel \"log(time/s)\"\n")
	fout.write("plot 'fort.98' u 1:2 w l t \""+specs[0].name+"\",\\\n")
	for i in range(neq-2):
		fout.write("'' u 1:"+str(i+3)+" w l t \""+specs[i+1].name+"\",\\\n")
	fout.write("'' u 1:"+str(neq+2)+" w l t \""+specs[neq-1].name+"\"\n")
	fout.close()

#*********************************
#*********************************
#*********************************
#copy other files to build
print "- copying others...",

#copy static files to build
if(is_test):
	print "- copying test to /build...",
	if(useDvodeF90):
		shutil.copyfile("tests/MakefileF90", buildFolder+"Makefile")
	elif(buildCompact):
		shutil.copyfile("tests/MakefileCompact", buildFolder+"Makefile")
	else:
		shutil.copyfile("tests/Makefile", buildFolder+"Makefile")

	test_file = "tests/test_"+test_name+".f90"
	plot_file = "tests/plots/plot_"+test_name+".gps"
	#chech if test file exists
	try:
		with open(test_file) as f: pass
	except IOError as e:
		print
		print "ERROR: Test file \""+test_file+"\" doesn't exist!"
		sys.exit()
	shutil.copyfile(test_file, buildFolder+"test.f90")
	if(has_plot): shutil.copyfile(plot_file, buildFolder+"plot.gps")
	print " done!"

#copy solver files to build
print "- copying solver to /build...",
if(useDvodeF90):
	shutil.copyfile("solver/dvode_f90_m.f90", buildFolder+"dvode_f90_m.f90")
else:
	shutil.copyfile("solver/opkdmain.f", buildFolder+"opkdmain.f")
	shutil.copyfile("solver/opkda1.f", buildFolder+"opkda1.f")
	shutil.copyfile("solver/opkda2.f", buildFolder+"opkda2.f")
print " done!"
print
print "You'll find the necessary files in "+buildFolder
print "Example call to the solver in "+buildFolder+"test.f90"
print "Example Makefile in "+buildFolder+"Makefile"

#check for large reaction set
if(len(reacts)>500 and not(use_implicit_RHS)):
	print
	print "WARNING: "+str(len(reacts))+" reactions found! Using implicit RHS (option -iRHS) could be more efficient"
	print "	 and also allows faster compilation."

print
if(not(is_test)):
	print "Call KROME in your code as:"
	print "    call krome(x(:), gas_density, gas_temperature, time_step)"
	print "where x(:) is a real*8 array of size "+str(len(specs)-4)+(" of the mass fractions" if useX else " of number densities [1/cm3]")
	print "gas_density " + ("is the gas density in [g/cm3]" if useX else "is a dummy argument")
	print "time_step is the integration time-step in [s]"
else:
	print "This is a test. To run it just type:"
	print "> cd build/"
	print "> make"
	print "> ./krome"
print
print "Everything done, goodbye!"


get_quote()

