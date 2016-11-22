import utils,itertools

class species():

	name = nameLatex = nameFile = nameHtml = nameHref = None

	#*****************
	#constructor
	def __init__(self,speciesName,atomSet):


		#check for upper case atoms, e.g. HE instead of He
		for atom in atomSet.keys():
			if((atom.upper()!=atom) and (atom.upper() in speciesName)):
				orgSpeciesName = speciesName
				speciesName = speciesName.replace(atom.upper(),atom)
				print "WARNING: "+atom+" found upper case in "+orgSpeciesName+"! Replaced as "+speciesName

		self.name = speciesName

		self.nameLatex = self.getLatexName()
		self.nameFile = self.nameHash = self.name.lower().replace("+","j").replace("-","w")
		self.nameHtml = self.getHtmlName()
		self.nameHref = self.getHrefName()

		me = 9.10938356e-28 #electron mass, g

		#store keys sorted by inverse length
		atoms = sorted(atomSet.keys(), key=lambda x:len(x),reverse=True)
		#produce unique character combinations
		alpha = ["".join(x) for x in list(itertools.product("XYZ",repeat=4))]
		#check to have enough combinations
		if(len(atoms)>len(alpha)): sys.exit("ERROR: in species parser alpha needs to be extended!")

		#local species name
		specName = speciesName #.upper()

		#compute charge
		positiveCount = specName.count("+")
		negativeCount = specName.count("-")
		self.charge = positiveCount-negativeCount
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
		self.explodedFull = self.exploded + (["+"]*positiveCount) + (["-"]*negativeCount)

		#store atoms
		nonAtoms = ["+","-"]
		self.atoms = [x for x in self.exploded if not(x in nonAtoms) and not(utils.isNumber(x))]

		#compute mass using dictionary as reference
		self.mass = sum([atomSet[x] for x in self.exploded])
		if(speciesName.upper()!="E"): self.mass -= self.charge*me

	#**********************
	def getHtmlName(self):
		name = list(self.name)
		if(self.name.upper()=="E"): return "e<sup>-</sup>"
		nameHtml = []
		for x in name:
			if(utils.isNumber(x)):
				y = "<sub>"+x+"</sub>"
			elif(x=="+" or x=="-"):
				y = "<sup>"+x+"</sup>"
			else:
				y = x
			nameHtml.append(y)
		return ("".join(nameHtml))

	#******************
	def getHrefName(self):
		url = "species_"+str(self.nameFile)+".html"
		return "<a href=\""+url+"\">"+self.getHtmlName()+"</a>"

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

	#*****************
	#get "engineered" enthalpy kJ/mol
	def getEnthalpy(self,thermochemicalData,Tgas=298.15):


		#gas constant kJ/mol/K
		Rgas = 8.3144598

		#use so-called electron convention
		if(self.name=="E"): return 5./2.*Rgas*Tgas

		#CRs has no enthalpy
		if(self.name=="CR"): return 0e0

		#return None if unkonw element
		if(not(self.name in thermochemicalData)): return None

		myData = thermochemicalData[self.name]
		if(Tgas>myData["Tmid"]):
			coefs = myData["highT"]
		else:
			coefs = myData["lowT"]

		HRT = coefs[0] + \
			coefs[1]*Tgas/2. + \
			coefs[2]*Tgas**2/3. + \
			coefs[3]*Tgas**3/4. + \
			coefs[4]*Tgas**4/5. + \
			coefs[5]/Tgas

		#kJ/mol
		return HRT*Rgas*Tgas/1e3

	#*********************
	def makeHtmlPage(self,myNetwork):

		fname = "htmls/species_"+str(self.nameFile)+".html"

		tableHeader = "<tr>"+("<th>"*30)

		tableFormation = []
		tableDestruction = []
		for reactions in myNetwork.reactions:
			if(self.name in [x.name for x in reactions.products]):
				tableFormation.append(reactions)
			if(self.name in [x.name for x in reactions.reactants]):
				tableDestruction.append(reactions)

		fout = open(fname,"w")
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">"+self.nameHtml+"</p>\n")
		fout.write("<br>\n")
		fout.write("<a href=\"indexSpecies.html\">back</a>\n")
		fout.write("<p>Enthalpy @ 298.15K: <b>" \
			+ str(self.getEnthalpy(myNetwork.thermochemicalData)) + "</b> kJ/mol</p>")

		fout.write("<br><br>\n")
		fout.write("<p style=\"font-size:20px\">Formation channels</p><br>\n")
		fout.write("<table width=\"50%\">\n")
		fout.write(tableHeader+"\n")
		icount = 0
		for reaction in tableFormation:
			bgcolor = ""
			if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
			fout.write("<tr bgcolor=\""+bgcolor+"\" valign=\"baseline\">"+reaction.getReactionHtmlRow(self)+"\n")
			icount += 1
		fout.write(tableHeader+"\n")
		fout.write("</table>\n")

		fout.write("<br><br>\n")
		fout.write("<p style=\"font-size:20px\">Destruction channels</p><br>\n")
		fout.write("<table width=\"50%\">\n")
		fout.write(tableHeader+"\n")
		icount = 0
		for reaction in tableDestruction:
			bgcolor = ""
			if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
			fout.write("<tr bgcolor=\""+bgcolor+"\" valign=\"baseline\">"+reaction.getReactionHtmlRow(self)+"\n")
			icount += 1
		fout.write(tableHeader+"\n")
		fout.write("</table>\n")

		fout.write(utils.getFooter("footer.php"))
		fout.close()





