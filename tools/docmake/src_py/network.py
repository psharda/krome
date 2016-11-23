import sys,glob,os
import reaction,utils,options

class network:

	species = None
	speciesDictionary = None

	#************************
	#class constructor
	def __init__(self,myOptions):

		#default format
		reactionFormat = "@format:idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"

		#get network file name
		fileName = myOptions.network

		#read atoms
		atomSet = utils.getAtomSet("atomlist.dat")
		#read shortcuts
		shortcuts = utils.getShortcuts()

		#load thermochemical data
		self.thermochemicalData = utils.getThermochemicalData("thermo30.dat")

		inBlockCR = False

		self.reactions = []
		print "reading network "+fileName
		#open file to read network
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue

			#get format if token found
			if(srow.startswith("@format:")):
				reactionFormat = srow
				continue

			#get extra variable definition if token found
			if(srow.startswith("@var:")):
				(variable,expression) = [x.strip() for x in srow.replace("@var:","").split("=")]
				shortcuts[variable] = expression

			if(srow.lower().startswith("@cr_start") or srow.lower().startswith("@cr_begin")): inBlockCR = True
			if(srow.lower().startswith("@cr_stop") or srow.lower().startswith("@cr_end")): inBlockCR = False

			#change reaction type to CR
			reactionType = "standard" #default
			if(inBlockCR): reactionType = "CR"

			#skip other tokens
			if(srow.startswith("@")): continue

			#parse row line for reaction
			myReaction = reaction.reaction(srow,reactionFormat,atomSet,reactionType)

			#add parsed reaction to reactions structure in network
			self.reactions.append(myReaction)
		fh.close()

		#merge reactions with multiple rates
		print "merging multiple reactions"
		self.mergeReactions()

		#get a list with missing reactions (reactants)
		self.missingReactions = self.getMissing()

		#get ranges from option file
		varRanges = dict()
		for rng in myOptions.range:
			#store range name and range limits
			(rangeName,rangeValue) = [x.strip() for x in rng.split("=")]
			varRanges[rangeName] = [float(x) for x in rangeValue.split(",")]
			print rangeName, varRanges[rangeName]

		#clear folders (html and png)
		self.clearFolders()

		#prepare graphs
		self.makeGraph()

		#create reaction index page
		self.makeHtmlIndex()
		self.makeHtmlMissingReactionIndex(myOptions)
		self.makeHtmlReactionIndex(myOptions)
		self.makeHtmlSpeciesIndex()
		self.makeHtmlGraphIndex()
		self.makeHtmlAllRates(myOptions)


		#prepare html pages for species
		for mySpecies in self.getSpecies():
			mySpecies.makeHtmlPage(self)

		#plotting rates
		print "plotting rates..."
		icount = 0
		#loop on reactions to evaluate and plot
		for myReaction in self.reactions:
			print str(int(icount*1e2/len(self.reactions)))+"%", myReaction.getVerbatim()
			myReaction.plotRate(shortcuts,varRanges)
			myReaction.makeHtmlPage(myOptions)
			icount += 1

	#**************
	#get all network species
	def getSpecies(self):
		if(self.species!=None): return self.species
		self.species = []
		#loop on reactions
		for reaction in self.reactions:
			self.species += reaction.reactants + reaction.products

		#unique list of species
		dicSpec = dict()
		for mySpecies in self.species:
			dicSpec[mySpecies.nameHash] = mySpecies
		self.species = dicSpec.values()

		return self.species

	#**************
	#get all network species in a dictionary with NAME as keys
	def getSpeciesDictionary(self):
		if(self.speciesDictionary!=None): return self.speciesDictionary

		#unique list of species
		self.speciesDictionary = dict()
		for mySpecies in self.getSpecies():
			self.speciesDictionary[mySpecies.name] = mySpecies

		return self.speciesDictionary

	#*******************
	def getAtoms(self):
		atoms = []
		for mySpecies in self.getSpecies():
			atoms += mySpecies.atoms
		return list(set(atoms))

	#**************
	#clean temporary folders
	def clearFolders(self):
		#files to be deleted (folder:extension)
		folders = {"pngs":"png", \
			"htmls":"html", \
			"evals":"json", \
			"epss":"eps", \
			"dots":"dot"}

		#loop on folders and extensions
		for (path,extension) in folders.iteritems():
			removePath = path+"/*."+extension
			print "removing "+removePath
			#get files list
			filelist = glob.glob(removePath)
			#remove files
			for fname in filelist:
				os.remove(fname)

	#**************
	#merge reactions with multiple rates
	def mergeReactions(self):
		reactionList = dict()
		#merge reactions
		for myReaction in self.reactions:
			myHash = myReaction.getReactionHash()
			if(myHash in reactionList):
				reactionList[myHash].merge(myReaction)
			else:
				reactionList[myHash] = myReaction
		self.reactions = [v for (k,v) in reactionList.iteritems()]


	#***************
	def makeGraph(self):

		import subprocess

		#loop on all available atoms in the network
		for refAtom in self.getAtoms():
			print "creating graphs for "+refAtom+"..."

			#connection dictionary
			networkMap = dict()
			#loop on reactions
			for myReaction in self.reactions:
				#loop on reactants
				for myR in myReaction.reactants:
					#skip species if reference atom is not contained
					if(not(refAtom in myR.atoms) and refAtom!=None): continue
					#get reactant partners in this reaction (other reactants)
					partners = [x.name for x in myReaction.reactants if (x.name!=myR.name)]
					#loop on products
					for myP in myReaction.products:
						if(myP.name==myR.name): continue
						#skip species if reference atom is not contained
						if(not(refAtom in myP.atoms) and refAtom!=None): continue
						#create dictionary if not present
						if(not(myR.name in networkMap)):
							networkMap[myR.name] = dict()
						#create list if not present
						if(not(myP.name in networkMap[myR.name])):
							networkMap[myR.name][myP.name] = []
						#add partners to the matrix
						networkMap[myR.name][myP.name] += partners

			#dot, png, eps files
			fname = "dots/network_atom_"+str(refAtom)+".dot"
			fnamePNG = "pngs/network_atom_"+str(refAtom)+".png"
			fnameEPS = "epss/network_atom_"+str(refAtom)+".eps"

			#write dot
			fout = open(fname,"w")
			fout.write("digraph G{\n")

			if(len(networkMap.keys())==0): continue

			nodeLabels = dict()
			#loop on starting edges
			for myR in networkMap.keys():
				nodeLabels[myR.lower()] = myR
				#loop on ending edges (and partners)
				for (myP,partners) in networkMap[myR].iteritems():
					nodeLabels[myP.lower()] = myP
					labelList = sorted(list(set(partners)))
					label = ",".join(labelList)
					#write edge to dot
					fout.write("\""+myR.lower()+"\" -> \""+myP.lower()+"\" [label=\""+label+"\"];\n")

			for (k,v) in nodeLabels.iteritems():
				fout.write("\""+k+"\" [label=\""+v+"\"];\n")

			fout.write("}\n")
			fout.close()

			#write PNG with dot
			myCall = "dot "+fname+" -Tpng -o "+fnamePNG
			subprocess.call(myCall.split(" "))

			#write EPS with dot
			myCall = "dot "+fname+" -Teps -o "+fnameEPS
			subprocess.call(myCall.split(" "))


	#****************
	#get a list of the missing reactions (reactants only)
	def getMissing(self):

		#get list of species
		species = self.getSpecies()

		#evaluate expected reactions
		print "expected reactions:",len(species)**2+len(species)

		#unimolecular reactants to be skipped
		skipUnimol = ["CR","E"]

		allReactants = []
		#loop on species (reactant 1)
		for species1 in species:
			#add unimolecular (only neutral or anion)
			if((species1.charge<1) and not(species1.name in skipUnimol)): allReactants.append([species1.name])
			#loop on species (reactant 2) for bimolecular
			for species2 in species:
				#exclude repulsive
				if(species1.charge*species2.charge>0): continue

				#exclude double ionization CR
				hasCR = (species1.name=="CR" or species2.name=="CR")
				hasElectron = (species1.name=="E" or species2.name=="E")
				hasCation = (species1.charge>0 or species2.charge>0)
				if(hasCation and hasCR): continue
				if(hasCR and hasElectron): continue

				#add bimolecular
				reactants = [species1.name,species2.name]
				allReactants.append(sorted(reactants))


		#get species dictionary to return species objects instead of names
		speciesDictionary = self.getSpeciesDictionary()

		#create a dictionary of products using exploded as key
		# e.g. -_C_H_H_O for H2O + C-
		productDictionary = dict()
		#loop on all reactants (products are the same since species1+species2 generated)
		for reactants in allReactants:
			reactExploded = []
			#loop on reactants for exploded (not explodedFull contains charge signs)
			for RR in reactants:
				reactExploded += speciesDictionary[RR].explodedFull
			#concatenate exploded, e.g. -_C_H_H_O for H2O + C-
			reactExp = ("_".join(sorted(reactExploded)))
			#remove mutual sign neutralization
			reactExp = reactExp.replace("+_-_","")
			#store products in a list (init list if not present)
			if(not(reactExp in productDictionary)):
				productDictionary[reactExp] = []
			#store unique products
			if(not(reactants in productDictionary[reactExp])): productDictionary[reactExp].append(reactants)


		#write number of potential reactions found
		print "expected reactions (after criteria):",len(allReactants)

		#store reactants for each reaction
		networkReactants = []
		for myReaction in self.reactions:
			reactants = sorted([x.name for x in myReaction.reactants])
			networkReactants.append(reactants)

		#search for missing rates
		missing = []
		#loop on all possible reactants groups
		for reactant in allReactants:
			found = False
			#loop on network reactant groups
			for rectantNtw in networkReactants:
				if(reactant==rectantNtw):
					found = True
					break
			#store not found if not present already
			if(not(found) and not(reactant in missing)): missing.append(reactant)

		#sort missing reactions by reactants name
		missing = sorted(missing,key=lambda x:("_".join(x)))
		print "missing reactions:",len(missing)

		missingReactions = []
		#loop on missing reactions to search for possible branches
		for reactants in missing:
			#create exploded reaction
			reactExploded = []
			for RR in reactants:
				reactExploded += speciesDictionary[RR].explodedFull
			reactExp = ("_".join(sorted(reactExploded)))
			#remove mutual sign neutralization
			reactExp = reactExp.replace("+_-_","")
			#skip reactants==products, e.g. H+OH -> H+OH (and convert names into species objects)
			productsList = [[speciesDictionary[x] for x in products] for products in productDictionary[reactExp] if(products!=reactants)]
			#convert names in to species objects
			reatantsList = [speciesDictionary[x] for x in reactants]
			missingReactions.append({"reactants":reatantsList,"products":productsList})

		#return list of objects reactants
		return missingReactions


	#****************
	def makeHtmlIndex(self):
		fname = "htmls/index.html"

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">KROME main menu</p>\n")
		fout.write("<br><br>")
		fout.write("<ul>")
		fout.write("<li><a href=\"indexReactions.html\">Reactions</a></li>")
		fout.write("<li><a href=\"indexSpecies.html\">Species</a></li>")
		fout.write("<li><a href=\"indexGraph.html\">Graphs</a></li>")
		fout.write("<li><a href=\"indexMissingReactions.html\">Missing reactions</a></li>")
		fout.write("</ul>")
		fout.write(utils.getFooter("footer.php"))
		fout.close()

	#****************
	def makeHtmlGraphIndex(self):

		fname = "htmls/indexGraph.html"

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">Graphs</p>\n")
		fout.write("<a href=\"index.html\">back</a><br>\n")
		fout.write("<br><br>\n")

		#loop on all available atoms in the network
		for refAtom in self.getAtoms():

			fnamePNG = "pngs/network_atom_"+str(refAtom)+".png"
			fnameEPS = "epss/network_atom_"+str(refAtom)+".eps"
			if(not(os.path.isfile(fnamePNG))): continue

			fout.write("<p style=\"font-size:20px\">Graph for "+refAtom+"</p>\n")
			fout.write("<p>&nbsp;Download <a href=\"../"+fnamePNG+"\" target=\"_blank\">PNG</a>")
			fout.write(" <a href=\"../"+fnameEPS+"\" target=\"_blank\">EPS</a></p><br>\n")
			fout.write("<a target=\"_blank\" href=\"../"+fnamePNG+"\">")
			fout.write("<img src=\"../"+fnamePNG+"\" style=\"{max-width:500px;}\"></a>\n")
			fout.write("<br><br>\n")

		#add footer
		fout.write(utils.getFooter("footer.php"))
		fout.close()


	#****************
	#create reaction list index as html page
	def makeHtmlSpeciesIndex(self):

		fname = "htmls/indexSpecies.html"

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">Species</p>\n")
		fout.write("<a href=\"index.html\">back</a><br>\n")
		fout.write("<br><br>\n")
		#reaction table
		fout.write("<table>\n")
		fout.write("<tr><th><th>\n")
		fout.write("<tr><td>name<td>&Delta;H (K)\n")
		fout.write("<tr><th><th>\n")
		icount = 0
		#loop on reactions
		for mySpecies in sorted(self.getSpecies(),key=lambda x:x.name):
			bgcolor = ""
			if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
			enthalpy = mySpecies.getEnthalpy(self.thermochemicalData)
			fout.write("<tr bgcolor=\""+bgcolor+"\"><td>&nbsp;"+mySpecies.getHrefName()+"&nbsp;<td>"+str(enthalpy)+"\n")
			icount += 1
		fout.write("<tr><th><th>\n")
		fout.write("</table>\n")

		#add footer
		fout.write(utils.getFooter("footer.php"))
		fout.close()


	#****************
	#create reaction list index as html page
	def makeHtmlReactionIndex(self,myOptions):

		fname = "htmls/indexReactions.html"

		tableHeader = "<tr>"+("<th>"*30)

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">Reactions</p>\n")
		fout.write("<a href=\"index.html\">back</a><br>\n")

		for variable in myOptions.getRanges().keys():
			fout.write("<a href=\"allRates_"+variable+".html\">All rates with <b>"+variable+"</b></a><br>\n")
		fout.write("<br><br>\n")
		#reaction table
		fout.write("<table width=\"50%\">\n")
		fout.write(tableHeader+"\n")
		icount = 0
		#loop on reactions
		for myReaction in self.reactions:
			bgcolor = ""
			if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
			fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+myReaction.getReactionHtmlRow()+"\n")
			icount += 1
		fout.write(tableHeader+"\n")
		fout.write("</table>\n")

		#add footer
		fout.write(utils.getFooter("footer.php"))
		fout.close()

	#****************
	#create reaction list index as html page
	def makeHtmlMissingReactionIndex(self,myOptions):

		fname = "htmls/indexMissingReactions.html"

		kJmol2K = 120.274 #kJ/mol->K

		tableHeader = "<tr>"+("<th>"*30)

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">Missing reactions*</p>\n")
		fout.write("<p style=\"font-size:10px\">*if reactants are present it doesn't check for missing branches!</p>\n")
		fout.write("<a href=\"index.html\">back</a><br>\n")


		rtypes = {1:"unimolecular", 2:"bimolecular", 3:"3-body"}
		for (nreact,rname) in rtypes.iteritems():
			fout.write("<br><br>\n")
			fout.write("<p style=\"font-size:20px\">"+rname.title()+", &Delta;H/K</p>\n")
			#reaction table
			fout.write("<table width=\"60%\">\n")
			fout.write(tableHeader+"\n")
			icount = 0
			#loop on reactions
			for reaction in self.missingReactions:
				#get reactants
				reactants = reaction["reactants"]
				#consider only number of reactants for the current reaction type
				if(len(reactants)!=nreact): continue
				#background color
				bgcolor = utils.getHtmlProperty("tableRowBgcolor")
				#create html row
				row = "<td>"+str(icount+1)+"<td>"+("<td>+<td>".join([x.getHtmlName() for x in reactants]))
				#if no branches available (given the network, skip)
				if(len(reaction["products"])==0): continue

				reactantsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in reactants]

				#add an arrow at the end of the row
				row += "<td><td>&rarr;"
				#write row to file
				fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+row+"\n")
				#count <td> elements
				ntds = row.count("<td>")
				#loop on products
				for products in reaction["products"]:
					productsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in products]
					bgcolor = "" #this row has no bgcolor
					#create row from products name
					rowx = "<td>&rarr;<td>"+("<td>+<td>".join([x.getHtmlName() for x in products]))
					#add n-1 td as offset
					row = ("<td>"*(ntds-1))+rowx+("<td>"*(10-rowx.count("<td>")))
					if(None in (productsEnthalpy+reactantsEnthalpy)):
						row += "<td>missing enthalpy data"
					else:
						DeltaH = sum(productsEnthalpy)-sum(reactantsEnthalpy)
						DeltaH_K = kJmol2K*DeltaH
						row += "<td>"+utils.htmlExp(DeltaH_K)
						if(DeltaH_K<0e0):
							row += "<td>&#10004;&#10004;"
						elif(DeltaH_K>=0 and DeltaH_K<1e4):
							row += "<td>&#10004;"
						else:
							row += "<td>&#10006;"
					#write html row to file
					fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+row+"\n")
				icount += 1

			fout.write(tableHeader+"\n")
			fout.write("</table>\n")

		#add footer
		fout.write(utils.getFooter("footer.php"))
		fout.close()

	#*******************
	def makeHtmlAllRates(self,myOptions):

		#loop on variables
		for variable in myOptions.getRanges().keys():

			#prepare a file for each variable
			fname = "htmls/allRates_"+variable+".html"

			#open file to write
			fout = open(fname,"w")
			#add header
			fout.write(utils.getFile("header.php"))
			fout.write("<p style=\"font-size:30px\">All rate plots with <b>"+variable+"</b></p>\n")
			fout.write("<a href=\"indexReactions.html\">back</a><br>\n")
			fout.write("<br><br>\n")
			#reaction table
			fout.write("<table>\n")
			icount = 0
			#loop on reactions
			for myReaction in self.reactions:
				#skip when the variable is not in the reaction rate
				if(not(myReaction.hasVariable(myOptions,variable))): continue
				if(icount%1==0): fout.write("<tr><td>")
				fnamePNG = "../pngs/rate_"+str(myReaction.getReactionHash())+"_"+variable+".png"
				linkURL = "<a href=\"rate_"+myReaction.getReactionHash()+".html\">details</a>"
				fout.write("<img src=\""+fnamePNG+"\" alt=\"&#9888; MISSING: "+myReaction.getVerbatim()+"\"><br>"+linkURL)
				fout.write("<tr height=\"10px\"><td>")
				icount += 1
			fout.write("</table>\n")

			#add footer
			fout.write(utils.getFooter("footer.php"))
			fout.close()



