import utils,itertools,os,sys

class species():

	name = nameLatex = nameFile = nameHtml = nameHref = None

	#*****************
	#constructor
	def __init__(self,speciesName,atomSet):


		#check for upper case atoms, e.g. HE instead of He
		# only if species (without signs) is not in the list
		if(not(speciesName.replace("+","").replace("-","") in atomSet.keys())):
			#loop on "atoms" names
			for atom in atomSet.keys():
				#if found lowercase replace it
				if((atom.upper()!=atom) and (atom.upper() in speciesName)):
					orgSpeciesName = speciesName
					speciesName = speciesName.replace(atom.upper(),atom)
					print "WARNING: "+atom+" found upper case in "+orgSpeciesName+"! Replaced as "+speciesName


		if(speciesName.startswith("J")): speciesName = speciesName[1:]+"_dust"

		#replace linear and chain prefix with underscore
		isoreplace = {"c-":"c_", "l-":"l_"}
		for (iso,isorep) in isoreplace.iteritems():
			if(speciesName.startswith(iso) and len(speciesName)>len(iso)):
				speciesName = speciesName.replace(iso,isorep)


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
		try:
			self.exploded = [atoms[alpha.index(x)] for x in exploded]
		except:
			print "ERROR: wanted to parse ",self.name
			print " but something went wrong with ", exploded
			print " Available atoms are:", atoms
			print " Add to atomlist.dat if needed."
			sys.exit()
		self.explodedFull = self.exploded + (["+"]*positiveCount) + (["-"]*negativeCount)

		#store atoms
		nonAtoms = ["+","-"]
		self.atoms = [x for x in self.exploded if(not(x in nonAtoms) and not(utils.isNumber(x)))]

		#compute mass using dictionary as reference
		self.mass = sum([atomSet[x] for x in self.exploded])
		if(speciesName.upper()!="E"): self.mass -= self.charge*me

		#init photochem variables
		self.xsecs = dict()
		self.phrates = dict()

		#load xsec from file
		self.loadXsecLeiden()

		#compute photo rates
		self.computePhRates()


	#**********************
	#load xsecs, structure is a dict
	# database_name
	#   |---energy---[list of energies, eV]
	#   |---phd/phi---[list of xsecs, eV]
	def loadXsecLeiden(self):

		fname = "xsecs/"+self.name+".dat"
		if(not(os.path.exists(fname))): return

		self.xsecs["leiden"] = dict()

		clight = 2.99792458e10 #cm/s
		hplanck = 4.135667662e-15 #eV*s

		#loop on file rows
		for row in open(fname,"rb"):
			srow = row.strip()
			if(srow==""): continue
			#look for format
			if(srow.startswith("# wavelength")):
				fmt = [x for x in srow.split(" ") if(x!="")][1:]
				#init arrays
				for fm in fmt:
					self.xsecs["leiden"][fm] = []
			if(srow.startswith("#")): continue
			arow = [float(x) for x in srow.split(" ") if(x!="")]
			#loop on data
			for i in range(len(fmt)):
				self.xsecs["leiden"][fmt[i]].append(arow[i])

		self.xsecs["leiden"]["energy"] = [clight*hplanck/(wl*1e-7) for wl in self.xsecs["leiden"]["wavelength"]]

		#reverse data
		for (k,v) in self.xsecs["leiden"].iteritems():
			self.xsecs["leiden"][k] = v[::-1]

	#**********************
	#compute photo rates, structure is a dict
	# database_name
	#   |---phd/phi
	#          |---radiation_fields---(rate,1/s)
	def computePhRates(self):
		import numpy as np
		from math import pi
		from photo import Jdraine,JBB,intJdraine,intJBB

		hplanck = 4.135667662e-15 #eV*s

		#radiation types
		frads = {"Draine1978":Jdraine, \
			"BB@4e3K":JBB, \
			"BB@1e4K":JBB, \
			"BB@2e4K":JBB}

		#loop on databases
		for (db,data) in self.xsecs.iteritems():
			self.phrates[db] = dict()
			xdata = data["energy"]
			#loop on rates
			for (k,ydata) in data.iteritems():
				#skip these keys
				if(k in ["energy","wavelength","photoabsorption"]): continue
				if(len(ydata)==0): continue
				self.phrates[db][k] = dict()
				#loop on radiation fluxes
				for radName,Jfun in frads.iteritems():
					if(radName.startswith("BB@")):
						Tbb = float(radName.replace("BB@","").replace("K",""))
						fscale = intJdraine()/intJBB(Tbb)
					fdata = []
					#loop on energy range
					for i in range(len(xdata)):
						sigma = ydata[i] #cm2
						energy = xdata[i] #eV
						if(radName.startswith("BB@")):
							Jrad = fscale*Jfun(energy,Tbb) #eV/cm2/sr
						else:
							Jrad = Jfun(energy) #eV/cm2/sr
						fdata.append(Jrad*sigma/energy)

					#compute integral and rate, 1/s
					kph = 4e0*pi*np.trapz(fdata,xdata)/hplanck

					#store rate, 1/s
					self.phrates[db][k][radName] = kph


	#****************************
	#plot xsecs to PNG
	def plotXsec(self):

		pngFileName = "pngs/xsec_"+self.nameFile+".png"
		import matplotlib.pyplot as plt

		#turn off interactivity
		plt.ioff()

		#cancel current plot
		plt.clf()

		#set aesthetics
		plt.grid(b=True, color='0.65',linestyle='--')
		plt.xlabel("energy/eV")
		plt.ylabel("xsec/cm2")
		plt.title(self.name)

		hasPlot = False
		#loop on databases
		for (db,data) in self.xsecs.iteritems():
			xdata = data["energy"]
			#loop on processes
			for (k,ydata) in data.iteritems():
				#skip these keys
				if(k in ["energy","wavelength","photoabsorption"]): continue
				#if no data no need to plot
				if(len(ydata)>0):
					hasPlot = True
					plt.plot(xdata,ydata,"-",label=k+" ("+db+")")

		plt.legend(loc='best')
		if(hasPlot): plt.savefig(pngFileName, dpi=150)


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
	#get polarizability, cm3
	def getPolarizability(self,polarizabilityData):

		if(self.name in polarizabilityData): return polarizabilityData[self.name]
		return None

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
		import urllib

		fname = "htmls/species_"+str(self.nameFile)+".html"

		tableHeader = "<tr>"+("<th>"*30)

		#do xsec plot
		self.plotXsec()

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
		polarizability = self.getPolarizability(myNetwork.polarizabilityData)
		if(polarizability!=None): polarizability /= 1e-24 #cm3->AA3
		fout.write("<p>&alpha;: <b>" \
			+ str(polarizability) + "</b> &Aring;<sup>3</sup></p>")
		fnameURL = "species_allrates_"+str(self.nameFile)+".html"
		fout.write("<p><a href=\""+fnameURL+"\">All plots in a single page</a></p>")

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
		fout.write("</table><br><br>\n")

		#put xsec png if data exists
		if(len(self.xsecs)>0):
			xsecPNG = "../pngs/xsec_"+self.nameFile+".png"
			fout.write("<img src=\""+xsecPNG+"\" width=\"500px\">\n")

			rows = []
			for (db,data) in self.phrates.iteritems():
				#loop on rates
				for (k,ydata) in data.iteritems():
					for radName,kph in ydata.iteritems():
						tr = "<td>&nbsp;"+k+"<td>"+db+"<td>"+radName+"<td>"+utils.htmlExp(kph)+"\n"
						rows.append([k+"_"+db+"_"+radName, tr])


			#table with integrated photorates
			thead = "<tr>"+"<th>"*4
			fout.write("<br><br>\n")
			fout.write("<p style=\"font-size:20px\">Photochemical rates (1/s)</p><br>\n")
			fout.write("<table width=\"50%\">\n")
			fout.write(thead+"\n")
			fout.write("<tr><td>&nbsp;process<td>database<td>radiation<td>rate (1/s)\n")
			fout.write(thead+"\n")
			icount = 0
			for row in sorted(rows,key=lambda x:x[0]):
				bgcolor = ""
				if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
				fout.write("<tr bgcolor=\""+bgcolor+"\">"+row[1])
				icount += 1
			fout.write(thead+"\n")
			fout.write("</table><br><br>\n")
			link = "http://home.strw.leidenuniv.nl/~ewine/photo/index.php?file=display_species.php&species="+urllib.quote(self.name, safe='')
			fout.write("<a href=\""+link+"\" target=\"_blank\">search on Leiden database</a>\n")

		fout.write(utils.getFooter("footer.php"))
		fout.close()

	#****************************
	def makeAllRatesHtmlPage(self,myNetwork,myOptions):

		fname = "htmls/species_allrates_"+str(self.nameFile)+".html"

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

		reactions = {"formation":tableFormation, "destruction":tableDestruction}

		#loop on reactions block
		for (k,data) in reactions.iteritems():
			fout.write("<br><br>\n")
			fout.write("<p style=\"font-size:20px\">"+k.title()+" channels</p><br>\n")
			#loop on reactions
			for reaction in data:
				#loop on variables
				for variable in myOptions.getRanges().keys():
					fnamePNG = "../pngs/rate_"+str(reaction.getReactionHash())+"_"+variable+".png"
					if(reaction.hasVariable(myOptions,variable)):
						#fout.write("<img src=\""+fnamePNG+"\">")
						linkURL = "<a href=\"rate_"+reaction.getReactionHash()+".html\">details</a> for "+reaction.getVerbatimHtml()
						fout.write("<img src=\""+fnamePNG+"\" width=\"700px\" alt=\"&#9888; MISSING: "+reaction.getVerbatim()+"\"><br>"+linkURL+"<br><br>")


		fout.write(utils.getFooter("footer.php"))
		fout.close()



