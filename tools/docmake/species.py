import utils,itertools

class species():

	name = nameLatex = nameFile = nameHtml = nameHref = None

	#*****************
	#constructor
	def __init__(self,speciesName,atomSet):

		self.name = speciesName

		self.nameLatex = self.getLatexName()
		self.nameFile = self.name.lower().replace("+","j").replace("-","w")
		self.nameHtml = self.getHtmlName()
		self.nameHref = self.getHrefName()

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

	#*********************
	def makeHtmlPage(self,myNetwork):

		fname = "htmls/species_"+str(self.nameFile)+".html"

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

		fout.write("<br><br>\n")
		fout.write("<p style=\"font-size:20px\">Formation channels</p>\n")
		fout.write("<table>\n")
		fout.write("<tr><th><th><th>\n")
		for reaction in tableFormation:
			fout.write("<tr valign=\"baseline\">"+reaction.getReactionHtmlRow(self)+"\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("</table>\n")

		fout.write("<br><br>\n")
		fout.write("<p style=\"font-size:20px\">Destruction channels</p>\n")
		fout.write("<table>\n")
		fout.write("<tr><th><th><th>\n")
		for reaction in tableDestruction:
			fout.write("<tr valign=\"baseline\">"+reaction.getReactionHtmlRow(self)+"\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("</table>\n")

		fout.write(utils.getFile("footer.php"))
		fout.close()





