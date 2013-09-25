import os,glob,shutil
from kromelib import *

class krome():
	#set defaults
	solver_MF = 222
	force_rwork = useHeating = doReport = checkConserv = useFileIdx = buildCompact = useEquilibrium = False
	use_implicit_RHS = use_photons = useTabs = useDvodeF90 = useTopology = useFlux = skipDup = False
	useCoolingAtomic = useCoolingH2 = useCoolingH2GP98 = useCoolingHD = useCoolingZ = use_cooling = useCoolingDust = False
	useCoolingCompton = useH2opacity = useCoolingCIE = False
	useCoolingZC = useCoolingZCp = useCoolingZSi = useCoolingZSip = useCoolingZO = useCoolingZOp = useCoolingZFe = useCoolingZFep = False
	useReverse = useCustomCoe = useODEConstant = cleanBuild = usePlainIsotopes = useDust = use_thermo = False
	usePhIoniz = useHeatingCompress = useHeatingPhoto = useHeatingChem = useDecoupled = useCoolingdH = useHeatingdH = False
	pedanticMakefile = useFakeOpacity = False
	useX = has_plot = doIndent = useTlimits = useODEthermo = True
	useDustGrowth = useDustSputter = useDustH2 = useDustT = False
	typeGamma = "DEFAULT"
	test_name = "default"
	is_test = False
	TlimitOpLow = "GE"
	TlimitOpHigh = "LT"
	customCoeFunction = "[CUSTOM COE FUNCTION NOT SET!]"
	buildFolder = "build/"
	TminAuto = 1e99
	TmaxAuto = -1e99
	RTOL = 1e-4 #default relative tolerance
	ATOL = 1e-20 #default absolute tolerance
	dustArraySize = dustTypesSize = 0
	dustTypes = []
	specs = []
	reacts = []
	dummy = molec()
	implicit_arrays = totMetals = ""
	thermodata = dict() #thermochemistry data (nasa polynomials)
	
	######################################
	#select test name
	def select_test(self,argv):
		for arg in argv:
			if("test=" in arg):
				self.is_test = True
				test_name = (arg.strip().replace("-test=",""))
				print "Reading option -test (test="+test_name+")"
				if(test_name=="planet"):
					[argv.append(x) for x in ["-useN","-reverse"]]
					filename = "networks/react_planet"
				elif(test_name=="cloud"):
					[argv.append(x) for x in ["-useN","-iRHS"]]
					filename = "networks/react_cloud"
				elif(test_name=="slowmanifold"):
					[argv.append(x) for x in ["-useN"]]
					filename = "networks/react_SM"
				elif(test_name=="shock1Dcool"):
					[argv.append(x) for x in ["-cooling=H2,HD,Z,DH"]]
					filename = "networks/react_primordial"
				elif(test_name=="shock1D"):
					filename = "networks/react_primordial"
				elif(test_name=="shock1Dphoto"):
					[argv.append(x) for x in ["-usePhIoniz","-heating=PHOTO","-cooling=ATOMIC,H2,HD,Z","-useEquilibrium"]]
					filename = "networks/react_primordial_photo"
				elif(test_name=="shock1Dlarge"):
					[argv.append(x) for x in ["-iRHS"]]
					filename = "networks/react_WH2008"
				elif(test_name=="dust"):
					[argv.append(x) for x in ["-dust=10,C,Si","-useN","-dustOptions=GROWTH,SPUTTER"]]
					filename = "networks/react_primordial"
				elif(test_name=="compact"):
					[argv.append(x) for x in ["-compact"]]
					filename = "networks/react_primordial"
				elif(test_name=="map"):
					[argv.append(x) for x in ["-cooling=ATOMIC,HD,H2", "-heating=PHOTO","-usePhIoniz"]]
					filename = "networks/react_primordial_photoH2"
				elif(test_name=="collapse"):
					[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE", "-heating=COMPRESS,CHEM"]]
					[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL"]]
					filename = "networks/react_primordial"
				elif(test_name=="collapseZ"):
					[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE,CII,OI", "-heating=COMPRESS,CHEM"]]
					[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL","-ATOL=1d-40","-forceMF=21"]]
					filename = "networks/react_primordialZ"
				elif(test_name=="collapseUV"):
					[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE", "-heating=COMPRESS,CHEM"]]
					[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL"]]
					filename = "networks/react_primordial_UV"
				elif(test_name=="collapseDUST"):
					[argv.append(x) for x in ["-cooling=ATOMIC,H2,COMPTON,CIE,DUST", "-heating=COMPRESS,CHEM"]]
					[argv.append(x) for x in ["-useH2opacity","-useN","-gamma=FULL","-dust=1,C","-dustOptions=TDUST"]]
					filename = "networks/react_primordial"
				elif(test_name=="kasting"):
					[argv.append(x) for x in ["-useN"]]
					filename = "networks/react_planet2"
				elif(test_name=="reverse"):
					[argv.append(x) for x in ["-useN","-reverse"]]
					filename = "networks/react_NO"
				elif(test_name=="atmosphere"):
					[argv.append(x) for x in ["-useN"]]
					filename = "networks/react_kast80"
				else:
					print "ERROR: test \""+test_name+"\" not present!"
					sys.exit()
				self.filename = filename
				self.test_name = test_name
				break

	##########################################
	def argparsing(self,argv):
		#you can select only one -forceMF
		if(("-forceMF=222" in argv) and ("-forceMF=21" in argv)):
			die("ERROR: options -forceMF=222 and -forceMF=21 are mutually exclusive: choose one.")

		#get filename
		if(not(self.is_test)): self.filename = argv[1]
		#chech if reactions file exists
		if(not(os.path.isfile(self.filename))): die("ERROR: Reaction file \""+self.filename+"\" doesn't exist!")


		#set implicit RHS
		if("-iRHS" in argv):
			self.use_implicit_RHS = True
			self.solver_MF = 222
			print "Reading option -iRHS"
		#force MF=21
		if("-forceMF=21" in argv):
			self.solver_MF = 21
			print "Reading option -forceMF=21"
		#force MF=222
		if("-forceMF=222" in argv):
			self.solver_MF = 222
			print "Reading option -forceMF=222"
		#use numeric density instead of fractions as input
		if("-useN" in argv):
			self.useX = False
			print "Reading option -useN"
		#force to use photons
		if("-usePhot" in argv):
			self.use_photons = True
			print "Reading option -usePhot"
		#use rate tables
		if("-useTabs" in argv):
			self.useTabs = True
			print "Reading option -useTabs"
		#use f90 solver
		if("-useDvodeF90" in argv):
			self.useDvodeF90 = True
			print "Reading option -useDvodeF90"
		#do report
		if("-report" in argv):
			self.doReport = True
			print "Reading option -report"
		#check mass conservation
		if("-checkConserv" in argv):
			self.checkConserv = True
			print "Reading option -checkConserv"
		#use reaction indexes in reaction file
		if("-useFileIdx" in argv):
			self.useFileIdx = True
			print "Reading option -useFileIdx"
		#write a single compact file krome_all.f90
		if("-compact" in argv):
			self.buildCompact = True
			print "Reading option -compact"
		#perform a clean build
		if("-clean" in argv):
			self.cleanBuild = True
			print "Reading option -clean"
		#build isotopes automatically
		if("-usePlainIsotopes" in argv):
			self.usePlainIsotopes = True
			print "Reading option -usePlainIsotopes"
		#use photoionization from Verner et al. 1996
		if("-usePhIoniz" in argv):
			self.usePhIoniz = True
			print "Reading option -usePhIoniz"
		#use equilibrium check to break loops earlier
		if("-useEquilibrium" in argv):
			self.useEquilibrium = True
			print "Reading option -useEquilibrium"
		#do not use temperature limits
		if("-noTlimits" in argv):
			self.useTlimits = False
			print "Reading option -noTlimits"
		#skip duplicated reactions
		if("-skipDup" in argv):
			self.skipDup = True
			print "Reading option -skipDup"
		#skip duplicated reactions
		if("-pedantic" in argv):
			self.pedanticMakefile = True
			print "Reading option -pedantic"
		#use reverse kinetics
		if("-reverse" in argv):
			self.useReverse = True
			print "Reading option -reverse"
		#use H2opacity following 
		if("-useH2opacity" in argv):
			self.useH2opacity = True
			print "Reading option -useH2opacity"

		#use cooling dT/dt in the ODE fex
		if("-skipODEthermo" in argv):
			self.useODEthermo = False
			print "Reading option -skipODEthermo"

		#determine Tgas limit operators
		for arg in argv:
			if("Tlimit=" in arg):
				self.myTlimit = (arg.strip().replace("-Tlimit=",""))
				self.myTlimit = myTlimit.replace("[","").replace("]","").split(",")
				self.TlimitOpHigh = myTlimit[1].strip().upper()
				self.TlimitOpLow = myTlimit[0].strip().upper()
				allOps = ["LE","LT","GE","GT"]
				if(not(self.TlimitOpLow in allOps) or not(self.TlimitOpHigh in allOps)):
					die("ERROR: on -Tlimit operators must be one of the followings: "+(", ".join(allOps)))
				print "Reading option -Tlimit (Low="+self.TlimitOpLow+", High="+self.TlimitOpHigh+")"
				break
		#determine cooling types
		for arg in argv:
			if("cooling=" in arg):
				myCools = arg.strip().replace("-cooling=","").split(",")
				allCools = ["ATOMIC","H2","HD","Z","DH","DUST","H2GP98","COMPTON","CIE",
						"CI","CII","SiI","SiII","OI","OII","FeI","FeII"]
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
				if("CIE" in myCools): self.useCoolingCIE = True
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

				self.use_thermo = True

				print "Reading option -cooling ("+(",".join(myCools))+")"
		
		#determine heating types
		for arg in argv:
			if("heating=" in arg):
				myHeat = arg.strip().replace("-heating=","").upper().split(",")
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
		for arg in argv:
			if("testFile=" in arg):
				myrwork = (arg.strip().replace("-testFile=",""))
				self.is_test = True
				self.test_name = myrwork
				print "Reading option -testFile (file="+str(myrwork)+")"
				break

		#force rwork size
		for arg in argv:
			if("forceRWORK=" in arg):
				myrwork = (arg.strip().replace("-forceRWORK=",""))
				self.force_rwork = True
				print "Reading option -forceRWORK (RWORK="+str(myrwork)+")"
				break
		#use custom function for coefficient instead of coe_tab()
		for arg in argv:
			if("useCustomCoe=" in arg):
				self.customCoeFunction = (arg.strip().replace("-useCustomCoe=",""))
				self.useCustomCoe = True
				print "Reading option -useCustomCoe (Expression="+str(customCoeFunction)+")"
				break
		#use function to append after each ODE
		for arg in argv:
			if("useODEConstant=" in arg):
				self.ODEConstant = (arg.strip().replace("-useODEConstant=",""))
				self.useODEConstant = True
				print "Reading option -useODEConstant (Constant="+str(ODEConstant)+")"
				break
		#dust
		hasDustOptions = False
		for arg in argv:
			if("dustOptions=" in arg): hasDustOptions = True
		for arg in argv:
			if("dust=" in arg):
				dustopt = (arg.strip().replace("-dust=",""))
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
				break
		#dust options
		for arg in argv:
			if("dustOptions=" in arg):
				if(not(self.useDust)): die("ERROR: you need -dust=[see help] to activate dust options!")
				dustopt = (arg.strip().replace("-dustOptions=",""))
				dustOptions = dustopt.split(",")
				if("GROWTH" in dustOptions): self.useDustGrowth = True
				if("SPUTTER" in dustOptions): self.useDustSputter = True
				if("H2" in dustOptions): self.useDustH2 = True
				if("TDUST" in dustOptions): self.useDustT = True
				print "Reading option -dustOptions (options="+(",".join(dustOptions))+")"
				break

		#project name folder
		for arg in argv:
			if("project=" in arg):
				self.projectName = (arg.strip().replace("-project=",""))
				print "Reading option -project (name="+str(projectName)+")"
				self.buildFolder = "build_"+projectName+"/"
				fout = open(projectName+".kpj","w")
				fout.write((" ".join(argv)))
				fout.close()
				break
		#typeGamma
		for arg in argv:
			if("gamma=" in arg):
				typeGamma = (arg.strip().replace("-gamma=",""))
				self.typeGamma = typeGamma.replace("\"","")
				print "Reading option -gamma (gamma="+str(self.typeGamma)+")"
				break

		#ATOL
		for arg in argv:
			if("ATOL=" in arg):
				self.ATOL = (arg.strip().replace("-ATOL=",""))
				print "Reading option -ATOL (ATOL="+str(self.ATOL)+")"
				break

		#RTOL
		for arg in argv:
			if("RTOL=" in arg):
				self.RTOL = (arg.strip().replace("-RTOL=",""))
				print "Reading option -RTOL (RTOL="+str(self.RTOL)+")"
				break

		#show help
		if("-help" in argv or "-h" in argv or "--help" in argv):
			get_usage()

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
				spec = arow[0].strip() #read species name
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
		hasFormat = False
		format_items = 4+len(ireact)+len(iprod)
		if(skipDup): fdup = open("duplicates.log","w")
		idxFound = tminFound = tmaxFound = rateFound = True
		specs = []
		reacts = []
		fh = open(filename,"rb") #OPEN FILE
		for row in fh:
			srow = row.strip() #stripped row
			if(srow.strip()==""): continue #looks for blank line
			if(srow[0]=="#"): continue #looks for comment line
			#search for format string
			if("@format:" in srow):
				idxFound = tminFound = tmaxFound = rateFound = False
				hasFormat = True
				srow = srow.replace("@format:","") #remove 
				print "Using custom format:"
				print srow
				arow = srow.split(",") #split format line
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
				if(not(rateFound)):
					print "ERROR: format must contain rate token"
					sys.exit()
			
				continue #skip format line (it is not a reaction line)
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
	###########################################add all the metals to cooling
	def addMetals(self):
		Zcools = []
		specs = self.specs
		if(self.useCoolingZC): Zcools.append("C") 
		if(self.useCoolingZCp): Zcools.append("C+") 
		if(self.useCoolingZO): Zcools.append("O") 
		if(self.useCoolingZOp): Zcools.append("O+") 
		if(self.useCoolingZSi): Zcools.append("Si") 
		if(self.useCoolingZSip): Zcools.append("Si+") 
		if(self.useCoolingZFe): Zcools.append("Fe") 
		if(self.useCoolingZFep): Zcools.append("Fe+")
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

	######################################
	def addDustT(self):
		specs = self.specs		
		dustTypes = self.dustTypes
		dustArraySize = self.dustArraySize
		dustTypesSize = self.dustTypesSize
		#add dust Temperatures as species
		if(self.useDustT):
			for dType in dustTypes:
				for i in range(dustArraySize):
						#create the object named dust_type_Tdust_index
						mymol = molec()
						mymol.name = "dust_"+dType+"_Tdust_"+str(i)
						mymol.charge = 0
						mymol.mass = 0.e0
						mymol.ename = "dust_"+dType+"_Tdust_"+str(i)
						mymol.fidx = "idx_dust_"+dType+"_Tdust_"+str(i)
						mymol.idx = len(specs)+1
						specs.append(mymol)
			print "Dust Tempperature added:",dustArraySize*dustTypesSize
		self.specs = specs

	###################################
	def addSpecial(self):
		specs = self.specs		
		#look for photons and CR to add to the species list
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
		#evaluate the number of chemical species (excluding, dust, dummies, Tgas, Tdust)
		self.nmols = len(self.specs)-4-self.dustArraySize*self.dustTypesSize
		self.ndustT = 0 #default for Tdust
		if(self.useDustT): 
			self.nmols -= self.dustArraySize*self.dustTypesSize #remove Tdust elements
			self.ndustT = self.dustArraySize*self.dustTypesSize #number of Tdust elements
		#print found values
		print
		print "ODEs needed:", len(self.specs)
		print "Reactions found:", len(self.reacts)
		print "Species found:", self.nmols


	###############################################
	def dumpNetwork(self):
		#dump species to log file
		fout = open("species.log","w")
		idx = 0
		for mol in self.specs:
			idx += 1
			fout.write(str(idx)+"\t"+mol.name+"\t"+mol.fidx+"\n")
		fout.close()

		#dump reactions to log file
		fout = open("reactions.log","w")
		idx = 0
		for rea in self.reacts:
			idx += 1
			fout.write(str(rea.idx)+"\t"+rea.verbatim+"\n")
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
				if(self.useDustT): sdustT += " + ["+str(dustArraySize)+" "+dType+"-dust T]"
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
						dns[nmols+j-1] += " - krome_dust_sput(Tgas,krome_dust_asize("
						dns[nmols+j-1] += str(j)+"),ntot,n("+str(nmols+j)+"))"
					if(self.useDustT):
						diffTdust = "( - n("+str(nmols+j)+") * krome_dust_asize3("+str(j)+") * krome_grain_rho &\n"
						diffTdust += " * dustEm(krome_dust_asize("+str(j)+"), n("+str(nmols+myndust+j)+")," 
						diffTdust += "dust_opt_asize_"+dType+", dust_opt_nu_"+dType+", dust_opt_Qabs_"+dType+") &\n"
						diffTdust += " + dustCool(krome_dust_asize2("+str(j)+"),n("+str(nmols+j)
						diffTdust += "), Tgas, n("+str(nmols+myndust+j)+"),ntot)) / boltzmann_erg / n("+str(nmols+j)+")"
						dns[nmols+myndust+j-1] += " + " + diffTdust

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
		for dn in dns:
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
				if((i+1) % 4 == 0): dns += "&\n" #break long lines
				if("-" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("-"," -"+str(RHSc[i])+".d0*")
				if("+" in RHSs[i] and RHSc[i]>1): dns += RHSs[i].replace("+"," +"+str(RHSc[i])+".d0*")
				if(RHSc[i]==1): dns+= " "+RHSs[i]
			dnw.append(dns)
			idn += 1
		self.dnw = dnw


	###############################################
	def createJAC(self):
		reacts = self.reacts
		specs = self.specs
		Tgas_species = self.Tgas_species
		ndust = self.dustArraySize * self.dustTypesSize
		#create explicit JAC for chemical species
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

		#create approximated Jacobian term for Tgas
		for i in range(neq):
			if(not(self.use_thermo)): break
			jsparse[i][Tgas_species.idx-1] = jsparse[Tgas_species.idx-1][i] = 1
			s = str(i+1)
			if(specs[i].name!="dummy" and specs[i].name!="CR" and specs[i].name!="g"):
				jac[Tgas_species.idx-1][i] = "if(abs(n("+s+") - jac_nold(" + s 
				jac[Tgas_species.idx-1][i] += "))>1d-10) pdj(idx_Tgas) = (jac_dn(idx_Tgas) - jac_dnold(idx_Tgas)) / (n("
				jac[Tgas_species.idx-1][i] += s + ") - jac_nold(" + s + "))"

		#TODO:create approximated Jacobian term for Tdust
		for i in range(neq):
			if(not(self.useDustT)): break
			for j in range(ndust):
				jTdust = nmols + (ndust) + j
				print i,j,jTdust
				sjTdust = str(jTdust)
				jsparse[i][jTdust] = jsparse[jTdust][i] = 1
				s = str(i+1)
				if(specs[i].name!="dummy" and specs[i].name!="CR" and specs[i].name!="g"):
					jac[jTdust][i] = "if(abs(n(" + s + ") - jac_nold(" + s
					jac[jTdust][i] += "))>1d-10) pdj(idx_Tgas) = (jac_dn(idx_Tgas) - jac_dnold(idx_Tgas)) / (n("
					jac[jTdust][i] += s + ") - jac_nold(" + s + "))"
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
		print


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

		print "LWM:",lwm,"LRW:",lrw
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
					fout.write("\tinteger,parameter::ntdust=" + str(ndust) + "\n")
					fout.write("\tinteger,parameter::ndustTypes=" + str(self.dustTypesSize) + "\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header())
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
		#********* CONSTANTS ****************
		fh = open("src/krome_constants.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_constants.f90","w")

		for row in fh:
			if(row[0]!="#"): fout.write(row)

		if(not(self.buildCompact)):
			fout.close()


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
				row = row.replace("#KROME_header", get_licence_header())

				if(row[0]!="#"): fouta.write(row)

			fouta.close()
			print "done!"
		else:
			print "WARNING: krome_user_commons.f90 not replaced!"

	###################################################
	def makeSubs(self):
		buildFolder = self.buildFolder
		reacts = self.reacts
		specs = self.specs
		thermodata = self.thermodata
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
					kstr = "\t" + sTlimit + " k("+str(x.idx)+") = " + x.krate + " !" + x.verbatim
					kstr = truncF90(kstr, 60,"*")
					fout.write(truncF90(kstr, 60,"/")+"\n")
			elif(srow == "#KROME_implicit_arrays"):
				fout.write(truncF90(self.implicit_arrays,60,","))
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
			elif(srow == "#KROME_var_reverse"):
				slen = str(len(specs))
				fout.write("real*8::p1("+slen+",7), p2("+slen+",7), Tlim("+slen+",3), p(7)\n")
			elif(srow == "#KROME_kc_reverse"):
				datarev = ""
				sp1 = sp2 = spt = ""
				print
				for x in specs:
					if(min(x.poly1)==0 and max(x.poly1)==0): continue
					sp1 += "p1("+x.fidx+",:)  = (/" + (", ".join([format_double(pp) for pp in x.poly1])) + "/)\n"
					sp2 += "p2("+x.fidx+",:)  = (/" + (", ".join([format_double(pp) for pp in x.poly2])) + "/)\n"
					spt += "Tlim("+x.fidx+",:)  = (/" + (", ".join([format_double(pp) for pp in x.Tpoly])) + "/)\n"
				fout.write(sp1+sp2+spt)
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header())
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
		#********* DUST ****************
		print "- writing krome_dust.f90...",
		fh = open("src/krome_dust.f90")
		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_dust.f90","w")


		dustPartnerIdx = dustQabs = dustOptInt = getTdust = ""
		if(self.useDust):
			itype = 0 #dust type index
			#loop in dust types
			for dType in self.dustTypes:
				itype += 1 #increase index
				dustPartnerIdx += "krome_dust_partner_idx("+str(itype)+") = idx_"+dType+"\n"
				if(self.useDustT):
					dustQabs += "call init_Qabs(\"opt"+dType+".dat\",dust_opt_Qabs_"+dType
					dustQabs += ",dust_opt_asize_"+dType+",dust_opt_nu_"+dType+")\n"
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
				fout.write(get_licence_header())
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
		if(self.useHeatingChem):
			Rref = []
			Pref = []
			kref = []
			Rref.append(sorted(["H","H","H"]))
			Pref.append(sorted(["H2","H"]))
			kref.append("4.48d0*n(idx_H)*n2H")
			Rref.append(sorted(["H2","H","H"]))
			Pref.append(sorted(["H2","H2"]))
			kref.append("4.48d0*n2H*n(idx_H2)")
			Rref.append(sorted(["H-","H"]))
			Pref.append(sorted(["H2","E"]))
			kref.append("3.53d0*n(idx_H)*n(idx_Hk)")
			Rref.append(sorted(["H2+","H"]))
			Pref.append(sorted(["H2","H+"]))
			kref.append("1.83d0*n(idx_H)*n(idx_H2j)")
			Rref.append(sorted(["H2","H"]))
			Pref.append(sorted(["H","H","H"]))
			kref.append("-4.48d0*n(idx_H)*n(idx_H2)")
			for rea in reacts:
				R = sorted([x.name for x in rea.reactants])
				P = sorted([x.name for x in rea.products])
				for i in range(len(Rref)):
					if(Rref[i]==R and Pref[i]==P):
						sclist = get_Tshortcut(rea,sclist) #get the shortcuts for temperature
						headchem = "!"+rea.verbatim + "\n"
						tklim = headchem + "if(Tgas." + rea.TminOp + "."  + rea.Tmin
						tklim += " .and. Tgas." + rea.TmaxOp + "." + rea.Tmax + ") then\n"
						HChem += tklim + "HChem = HChem + k("+str(rea.idx)+") * ("+kref[i] + ")\n end if\n\n"
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

				fout.write(get_licence_header())
			else:
				if(row.strip() == "#IFKROME_useHeatingdH" and (not(self.useHeatingdH) or len(dH_varsa)==0)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingCompress" and not(self.useHeatingCompress)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingPhoto" and not(self.useHeatingPhoto)): skip = True
				if(row.strip() == "#ENDIFKROME"): skip = False

				if(row.strip() == "#IFKROME_useHeatingChem" and not(self.useHeatingChem)): skip = True
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

		#*********ODE****************
		#write parameters in krome_ode.f90
		print "- writing krome_ode.f90...",

		fh = open("src/krome_ode.f90")

		if(self.buildCompact):
			fout = open(buildFolder+"krome_all.f90","a")
		else:
			fout = open(buildFolder+"krome_ode.f90","w")

		#root-finding based Tdust function (without dTdust/dt) PIPPO
		getTdust = ""
		#if(useDustT):
		#	iType = 0
		#	getTdust = "nd = "+str(dustArraySize)+" !number dust bins per types\n"
		#	for dType in dustTypes:
		#		getTdust += "call getTdust(krome_dust_T(nd*("+str(iType)+")+1:nd*"+str(iType+1)
		#		getTdust += "), dust_opt_Tbb_"+dType+", dust_opt_Em_"+dType+", dust_opt_asize_"+dType+",&\n"
		#		getTdust += " dust_opt_nu_"+dType+", dust_opt_Qabs_"+dType+", n(:))\n"
		#		getTdust = getTdust.replace("nd*(0)+1","1").replace("nd*1","nd").replace("nd*(1)","nd")
		#		iType += 1

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
						for x in dnw:
							fout.write("\t" + x + "\n")
			elif(srow == "#KROME_getTdust" and self.useDustT):
				fout.write(getTdust);
			elif(srow == "#KROME_dustSumVariables" and self.useDust):
				#add partner sum dust variable declarations
				dustSumVar = []
				for dType in dustTypes:
					dustSumVar.append("dSumDust"+dType)
				fout.write("\t real*8::" + (",".join(dustSumVar)) + "\n")
			elif(srow == "#KROME_header"):
				fout.write(get_licence_header())
			elif(srow == "#KROME_implicit_variables"):
				fout.write("real*8::rr\n")
				ris = (",".join(["r"+str(i+1) for i in range(self.maxnreag)]))
				pis = (",".join(["p"+str(i+1) for i in range(self.maxnprod)]))
				fout.write("integer::i,"+ris+","+pis+"\n")
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
							jacT = """!rough estimate for d^2T/dt/dx_i
								do i=1,neq-1
							         if(n(idx_Tgas)-jac_nold(idx_Tgas).ne.0.d0) then
								  pdj(i) = (jac_dnold(i) - jac_dn(i))/(jac_nold(idx_Tgas) - n(idx_Tgas))
								 end if
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
				fout.write(get_licence_header())
			elif(srow == "#KROME_common_alias"):
				fout.write("\tinteger,parameter::krome_nrea=" + str(len(reacts)) + "\n")
				fout.write("\tinteger,parameter::krome_nmols=" + str(nmols) + "\n")
				fout.write("\tinteger,parameter::krome_nspec=" + str(len(specs)) + "\n")
				fout.write("\tinteger,parameter::krome_ndust=" + str(dustArraySize*dustTypesSize) + "\n")
				if(self.useDustT): fout.write("\tinteger,parameter::krome_nTdust=" + str(dustArraySize*dustTypesSize) + "\n")
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
				fout.write("integer::"+(",".join(["r"+str(j+1) for j in range(self.maxnreag)]))+"\n")
			if(srow == "#KROME_arrs"):
				for j in range(self.maxnreag):
					fout.write("r"+str(j+1)+" = arr_r"+str(j+1)+"(i)\n")
			if(srow == "#KROME_arr_flux"):
				fout.write("arr_flux(i) = k(i)*"+("*".join(["n(r"+str(j+1)+")" for j in range(self.maxnreag)]))+"\n")

			if(row[0]!="#"): fout.write(row)
		if(not(self.buildCompact)):
			fout.close()

		print "done!"

	###############################
	def makeMain(self):
		buildFolder = self.buildFolder
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

			if(srow == "#IFKROME_useDustT" and not(self.useDustT)): skip = True
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
				fout.write(get_licence_header())
			elif(srow == "#KROME_rwork_array"):
				fout.write("\treal*8::rwork("+str(self.lrw)+")\n")
			elif(srow == "#KROME_iwork_array"):
				fout.write("\tinteger::iwork("+str(30+len(self.ia)+len(self.ja))+")\n")
			elif(srow == "#KROME_init_IAC"):
				fout.write("\t"+self.iaf+"\n")
			elif(srow == "#KROME_init_JAC"):
				fout.write("\t"+self.jaf+"\n")
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

		#copy static files to build
		if(self.is_test):
			print "- copying test to /build...",
			if(self.useDvodeF90):
				shutil.copyfile("tests/"+test_name+"/MakefileF90", buildFolder+"Makefile")
			elif(self.buildCompact):
				shutil.copyfile("tests/"+test_name+"/MakefileCompact", buildFolder+"Makefile")
			else:
				if(self.pedanticMakefile):
					shutil.copyfile("tests/Makefile_pedantic", buildFolder+"Makefile")
				else:
					shutil.copyfile("tests/"+test_name+"/Makefile", buildFolder+"Makefile")

			test_file = "tests/"+test_name+"/test.f90"
			plot_file = "tests/"+test_name+"/plot.gps"

			#chech if test file exists
			#if(not(os.path.isfile(filename))): die("ERROR: Test file \""+test_file+"\" doesn't exist!")

			shutil.copyfile(test_file, buildFolder+"test.f90")
			if(self.has_plot): shutil.copyfile(plot_file, buildFolder+"plot.gps")
			print " done!"

		#copy solver files to build
		print "- copying solver to /build...",
		if(self.useDvodeF90):
			shutil.copyfile("solver/dvode_f90_m.f90", buildFolder+"dvode_f90_m.f90")
		else:
			shutil.copyfile("solver/opkdmain.f", buildFolder+"opkdmain.f")
			shutil.copyfile("solver/opkda1.f", buildFolder+"opkda1.f")
			shutil.copyfile("solver/opkda2.f", buildFolder+"opkda2.f")
		print " done!"

		#kasting specific TODO:REPLACE
		if(test_name=="kasting"):
			shutil.copyfile("tests/kasting/init_specs", buildFolder+"init_specs")
			shutil.copyfile("tests/kasting/layers_data", buildFolder+"layers_data")


		#copy dust optical properties to build/
		if(self.useDustT):
			for dType in self.dustTypes:
				shutil.copyfile("data/opt"+dType+".dat", buildFolder+"opt"+dType+".dat")

	#######################################################
	def indent(self):
		buildFolder = self.buildFolder
		print "indenting...",
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
		if(not(self.is_test)):
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

		else:
			print "This is a test. To run it just type:"
			print "> cd build/"
			print "> make"
			print "> ./krome"
		print
		print "Everything done, goodbye!"
