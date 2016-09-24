import utils,itertools

class species():

	name = nameLatex = nameFile = nameHtml = None

	#*****************
	#constructor
	def __init__(self,speciesName,atomSet):

		self.name = speciesName

		self.nameLatex = self.getLatexName()
		self.nameHtml = self.getHtmlName()
		self.nameFile = self.name.replace("+","j").replace("-","w")

		me = 9.10938356e-28 #electron mass, g

		#store keys sorted by inverse length
		atoms = sorted(atomSet.keys(), key=lambda x:len(x),reverse=True)
		#produce unique character combinations
		alpha = ["".join(x) for x in list(itertools.product("XYZ",repeat=4))]
		#check to have enough combinations
		if(len(atoms)>len(alpha)): sys.exit("ERROR: in species parser alpha needs to be extended!")

		#species name uppercase
		specName = speciesName.upper()

		#compute charge
		self.charge = specName.count("+")-specName.count("-")
		if(speciesName.upper()=="E"): self.charge = -1
		#replace signs
		specName = specName.replace("+","").replace("-","")

		#replace atoms slash*separated
		for i in range(len(atoms)):
			specName = specName.replace(atoms[i],"/"+alpha[i]+"/")
		#replace double slashes
		while("//" in specName):
			specName = specName.replace("//","/")
		#split at slashes
		aspec = [x for x in specName.split("/") if x!=""]

		#search for number and when found multiplicate previous non-number
		exploded = []
		for a in aspec:
			if(utils.isNumber(a)):
				for j in range(int(a)-1):
					exploded.append(aold)
			else:
				exploded.append(a)
			aold = a

		#store exploded with real atom names
		self.exploded = [atoms[alpha.index(x)] for x in exploded]
		#compute mass using dictionary as reference
		self.mass = sum([atomSet[x] for x in self.exploded])
		if(speciesName.upper()!="E"): self.mass -= self.charge*me

	#**********************
	def getHtmlName(self):
		name = list(self.name)
		if(self.name.upper()=="E"): return "e<sup>-</sup>"
		nameLatex = []
		for x in name:
			if(utils.isNumber(x)):
				y = "<sub>"+x+"</sub>"
			elif(x=="+" or x=="-"):
				y = "<sup>"+x+"</sup>"
			else:
				y = x
			nameLatex.append(y)

		return ("".join(nameLatex))


	#**********************
	def getLatexName(self):
		name = list(self.name)
		if(self.name.upper()=="E"): return "e$^-$"
		nameLatex = []
		for x in name:
			if(utils.isNumber(x)):
				y = "$_"+x+"$"
			elif(x=="+" or x=="-"):
				y = "$^"+x+"$"
			else:
				y = x
			nameLatex.append(y)

		return ("".join(nameLatex))

