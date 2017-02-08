import sys,glob,os,shutil,hashlib
import reaction,utils,options

class network:

	species = None
	speciesDictionary = None

	#************************
	#class constructor
	def __init__(self,myOptions):


		#get network file name
		fileName = myOptions.network

		#read atoms
		atomSet = utils.getAtomSet("atomlist.dat")

		#load thermochemical data
		self.thermochemicalData = utils.getThermochemicalData("thermo30.dat")

		#read shortcuts
		shortcuts = utils.getShortcuts()

		#read network file in KROME format
		if(self.detectNetworkType(fileName)=="KIDA"):
			self.readFileKIDA(fileName,atomSet)
		else:
			self.readFileKROME(fileName,atomSet,shortcuts)

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
		self.backupEvaluationJSON()

		#prepare graphs
		self.makeGraph()


		#loop on reactions to evaluate
		for myReaction in self.reactions:
			myReaction.evaluateRate(shortcuts,varRanges)
			myReaction.makeHtmlPage(myOptions)

		#create reaction index page
		self.makeHtmlIndex()
		self.makeHtmlMissingBranchesIndex(myOptions)
		self.makeHtmlMissingReactionIndex(myOptions)
		self.makeHtmlReactionIndex(myOptions)
		self.makeHtmlSpeciesIndex()
		self.makeHtmlGraphIndex()
		self.makeHtmlAllRates(myOptions)
		self.makeHtmlMultipleRates(myOptions)

		self.deleteChangedPNGs(myOptions)

		#prepare html pages for species
		for mySpecies in self.getSpecies():
			mySpecies.makeHtmlPage(self)
			mySpecies.makeAllRatesHtmlPage(self,myOptions)

		#plotting rates
		print "plotting rates..."
		icount = 0
		#loop on reactions to plot
		for myReaction in self.reactions:
			print str(int(icount*1e2/len(self.reactions)))+"%", myReaction.getVerbatim()
			myReaction.plotRate(myOptions)
			icount += 1


	#********************
	def detectNetworkType(self,fileName):
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			if("," in srow): return "KROME"
		fh.close()
		return "KIDA"

	#************************
	def readFileKROME(self,fileName,atomSet,shortcuts):

		#default format
		reactionFormat = "@format:idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"

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



	#**************
	def readFileKIDA(self,fileName,atomSet):
		self.reactions = []
		print "reading network "+fileName

		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue

			reactionType = "KIDA"
			reactionFormat = ""

			#parse row line for reaction
			myReaction = reaction.reaction(srow,reactionFormat,atomSet,reactionType)

			#add parsed reaction to reactions structure in network
			self.reactions.append(myReaction)

		fh.close()

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
		#png and json are cleaned later if different
		folders = {"htmls":"html", \
			"epss":"eps", \
			"dots":"dot"}
		#"pngs":"png", \
		#"evals":"json", \


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
	#backup json files with rate evaluations (to avoid replotting)
	def backupEvaluationJSON(self):

		backupPath = "evals/*.json"
		#get files list
		filelist = glob.glob(backupPath)
		#copy files
		for fname in filelist:
			shutil.copyfile(fname, fname+".bak")

	#**************
	#delete PNG files with no json.bak file or different MD5
	# after this remove all json.bak files
	def deleteChangedPNGs(self,myOptions):

		#get file lists
		pngList = glob.glob("pngs/rate_*.png")

		#loop on variable ranges
		for rng in myOptions.range:
			#get range name
			rngName = rng.split("=")[0].strip()
			#produce end of the file name
			fileNameEnd = "_"+rngName+".png"
			#copy files
			for pngFname in pngList:
				survived = True
				#if png files end with range name
				if(fileNameEnd in pngFname):
					#create json and json.bak filenames
					jsonFname = pngFname.replace(fileNameEnd,"").replace("pngs/","evals/")+".json"
					jsonBakFname = jsonFname+".bak"
					#if bak file is missing remove png file
					if(not(os.path.exists(jsonBakFname))):
						#print "missing json bak", pngFname
						survived = False
						os.remove(pngFname)
					else:
						#check MD5 of json and json.bak
						md5 = hashlib.md5(open(jsonFname).read()).hexdigest()
						md5Bak = hashlib.md5(open(jsonBakFname).read()).hexdigest()
						#if md5 are different remove png file
						if(md5!=md5Bak):
							#print "diff MD5", pngFname,md5, md5Bak
							survived = False
							os.remove(pngFname)
					#if(survived): print "survived: "+pngFname

		#remove json.bak files
		bakList = glob.glob("evals/*.bak")
		for jsonBak in bakList:
			os.remove(jsonBak)

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

		#maximum number of nodes
		maxNodes = 50

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

			#skip network with too many nodes
			if(len(networkMap.keys())>maxNodes):
				print "WARNING: too many nodes, skipping graph!"
				continue

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
	#return all the possible reactants combinations with the current species
	# as an array, including unimol and bimol reactions
	def getAllReactantsCombinations(self,returnNames=True,includeRepulsive=False):

		#get list of species
		species = self.getSpecies()

		#unimolecular reactants to be skipped
		skipUnimol = ["CR","E"]

		allReactants = []
		#loop on species (reactant 1)
		for species1 in species:
			#add unimolecular (only neutral or anion)
			if((species1.charge<1) and not(species1.name in skipUnimol)):
				if(returnNames):
					allReactants.append([species1.name])
				else:
					allReactants.append([species1])
			#loop on species (reactant 2) for bimolecular
			for species2 in species:
				#exclude repulsive
				if(species1.charge*species2.charge>0 and not(includeRepulsive)): continue

				#exclude double ionization CR
				hasCR = (species1.name=="CR" or species2.name=="CR")
				hasElectron = (species1.name=="E" or species2.name=="E")
				hasCation = (species1.charge>0 or species2.charge>0)
				if(hasCation and hasCR): continue
				if(hasCR and hasElectron): continue

				#add bimolecular
				if(returnNames):
					reactants = [species1.name,species2.name]
					allReactants.append(sorted(reactants))
				else:
					allReactants.append([species1,species2])

		return allReactants


	#****************
	#get a list of the missing reactions (reactants only)
	def getMissing(self):

		#get list of species
		species = self.getSpecies()

		#limit to produce missing reactions
		expectedReactionsMax = int(1e4)

		#evaluate expected reactions
		expectedReactions = len(species)**2+len(species)
		print "expected reactions:", expectedReactions

		#skip if too many reactions
		if(expectedReactions>expectedReactionsMax):
			print "WARNING: too many expected reactions skipping missing reactions!"
			return []

		#get all possible reactants combinations
		allReactants = self.getAllReactantsCombinations()

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
			#create exploded reaction to compare and get branches
			reactExploded = []
			for RR in reactants:
				reactExploded += speciesDictionary[RR].explodedFull
			reactExp = ("_".join(sorted(reactExploded)))
			#remove mutual sign neutralization
			reactExp = reactExp.replace("+_-_","")
			#skip reactants==products, e.g. H+OH -> H+OH (and convert names into species objects)
			productsList = [[speciesDictionary[x] for x in products] for products in productDictionary[reactExp] if(products!=reactants)]
			productsListPrune = []
			for products in productsList:
				charges = [x.charge for x in products]
				#skip reaction that generates anion-cation products
				if((max(charges)+min(charges)==0) and (max(charges)>0)): continue
				productsListPrune.append(products)

			#convert names in to species objects
			reatantsList = [speciesDictionary[x] for x in reactants]
			missingReactions.append({"reactants":reatantsList,"products":productsListPrune})

		#return list of objects reactants
		return missingReactions

	#****************
	#given a list of species objects return the sum of their exploded
	# taking care of charges
	def getExplodedSpecies(self,speciesList):

		rexp = []
		#sum exploded
		for species in speciesList:
			rexp += species.explodedFull
		#wrap into underscores
		rexpFull = "_"+("_".join(sorted(rexp)))+"_"
		#remove +/- one by one until both are present
		while(("+" in rexpFull) and ("-" in rexpFull)):
			rexpFull = rexpFull.replace("+","",1)
			rexpFull = rexpFull.replace("-","",1)
		#remove +/E one by one until both are present
		while(("_+_" in rexpFull) and ("_E_" in rexpFull)):
			rexpFull = rexpFull.replace("+","",1)
			rexpFull = rexpFull.replace("E","",1)

		#for this hash E is just represented  as -
		rexpFull = rexpFull.replace("_E_","_-_")
		#remove CRs
		rexpFull = rexpFull.replace("CR","")
		#remove double underscores
		rexpFull = rexpFull.replace("__","_")

		#sort again replaced
		rexpFull = ("_".join(sorted(rexpFull.split("_"))))

		#remove leading and trailing underscores
		while(rexpFull.startswith("_")):
			rexpFull = rexpFull[1:]
		while(rexpFull.endswith("_")):
			rexpFull = rexpFull[:-1]

		return rexpFull


	#****************
	def getMissingBranch(self):

		missingBranch = []

		species = self.getSpecies()

		#get all possible species combinations
		allSpecies = dict()
		for spec1 in species:
			unimol = [spec1]
			rHash = self.getExplodedSpecies(unimol)
			if(not(rHash in allSpecies)): allSpecies[rHash] = []
			allSpecies[rHash].append(unimol)
			#bimolecular
			for spec2 in species:
				bimol = [spec1,spec2]
				rHash = self.getExplodedSpecies(bimol)
				if(not(rHash in allSpecies)): allSpecies[rHash] = []
				allSpecies[rHash].append(bimol)
				#trimolecular
				#for spec3 in species:
				#	trimol = [spec1,spec2,spec3]
				#	rHash = self.getExplodedSpecies(trimol)
				#	if(not(rHash in allSpecies)): allSpecies[rHash] = []
				#	allSpecies[rHash].append(trimol)

		#search all branches of a given reactions and store
		allReact = dict()
		for reaction in self.reactions:
			rHash = self.getExplodedSpecies(reaction.reactants)
			if(not(rHash in allReact)): allReact[rHash] = []
			allReact[rHash].append(reaction)


		#loop reactions to find missing and present branches
		for (rHash,reactions) in allReact.iteritems():

			#init data branch structure
			dataBranch = dict()
			dataBranch["reactants"] = reactions[0].reactants
			dataBranch["missingBranches"] = []
			dataBranch["presentBranches"] = []

			#store all reaction branches
			branches = []
			for reaction in reactions:
				branches.append(sorted([x.name for x in reaction.products]))

			#get reactant names (same for all reactions)
			reactNames = sorted([x.name for x in reactions[0].reactants])
			reactNamesNoCR = sorted([x.name for x in reactions[0].reactants if(x.name!="CR")])

			#get reaction hash for reactants only
			rHash = self.getExplodedSpecies(reactions[0].reactants)
			branchFound = []
			#loop on possible species combination for the given reactant hash
			for species in allSpecies[rHash]:
				branch = sorted([x.name for x in species])
				totalCharge = sum([x.charge for x in species])
				maxCharge = max([x.charge for x in species])
				hasElectron = ("E" in [x.name for x in species])
				#skip repulsive products
				if(totalCharge==0 and maxCharge!=0 and not(hasElectron)): continue
				#skip CR in products
				if("CR" in branch): continue
				#skip X+CR->X
				if(branch==reactNamesNoCR): continue
				#skip already-found branches
				if(branch in branchFound): continue
				#skip same reactant and products
				if(branch==reactNames): continue
				#store branches
				if(branch in branches):
					dataBranch["presentBranches"].append(species)
				else:
					dataBranch["missingBranches"].append(species)
				branchFound.append(branch)

			#copy data to structure
			missingBranch.append(dataBranch)

		return missingBranch


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
		fout.write("<li><a href=\"indexMissingBranches.html\">Missing branches</a></li>")
		fout.write("<li><a href=\"multipleRates.html\">Reactions with multiple rates</a></li>")
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
			fout.write("<img src=\"../"+fnamePNG+"\" width=\"700px\" style=\"{max-width:500px;}\"></a>\n")
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
		fout.write("<table width=\"60%\">\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("<tr><td>name<td>&Delta;H@0K (kJ/mol)<td>&Delta;H@298.15K (kJ/mol)\n")
		fout.write("<tr><th><th><th>\n")
		icount = 0
		#loop on reactions
		for mySpecies in sorted(self.getSpecies(),key=lambda x:x.name):
			bgcolor = ""
			if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
			enthalpy0 = mySpecies.getEnthalpy(self.thermochemicalData,Tgas=1e-40)
			enthalpy298 = mySpecies.getEnthalpy(self.thermochemicalData)
			fout.write("<tr bgcolor=\""+bgcolor+"\"><td>&nbsp;"+mySpecies.getHrefName()+"<td>"+str(enthalpy0)+"<td>"+str(enthalpy298)+"\n")
			icount += 1
		fout.write("<tr><th><th><th>\n")
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
		#reactions sorted by unsorted reaction hash (in the html list)
		reactionsSorted = sorted(self.reactions, key=lambda x:x.getReactionHashUnsorted())
		#loop on reactions
		for myReaction in reactionsSorted:
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
	def makeHtmlMissingBranchesIndex(self,myOptions):

		fname = "htmls/indexMissingBranches.html"

		kJmol2K = 120.274 #kJ/mol->K

		tableHeader = "<tr>"+("<th>"*30)

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">Missing branches, &Delta;H/K</p>\n")
		#fout.write("<p style=\"font-size:10px\">*if reactants are present it doesn't check for missing branches!</p>\n")
		fout.write("<a href=\"index.html\">back</a><br><br><br>\n")

		#sort reactions by the name of the first reactant (reactants are also sorted)
		sortedMissingBranch = sorted(self.getMissingBranch(),key=lambda x:sorted([sp.name for sp in x["reactants"]]))

		fout.write("<table width=\"60%\">\n")
		icount = 0
		for dataBranch in sortedMissingBranch:
			reactants = dataBranch["reactants"]
			reactantsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in reactants]
			row = "<td>"+str(icount+1)+"<td>"+("<td>+<td>".join(sorted([x.getHtmlName() for x in reactants])))

			row += ("<td>"*(20-row.count("<td>")))

			#background color
			bgcolor = utils.getHtmlProperty("tableRowBgcolor")

			#add an arrow at the end of the row
			row += "<td><td>&rarr;"
			#write row to file
			fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+row+"\n")

			#count <td> elements
			ntds = row.count("<td>")

			for listName in ["presentBranches","missingBranches"]:
				for products in dataBranch[listName]:
					bgcolor = "" #this row has no bgcolor
					productsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in products]
					bgcolor = "" #this row has no bgcolor
					#create row from products name
					rowx = "<td>&rarr;<td>"+("<td>+<td>".join(sorted([x.getHtmlName() for x in products])))
					#add n-1 td as offset
					row = ("<td>"*(ntds-1))+rowx+("<td>"*(10-rowx.count("<td>")))
					status = listName.replace("Branches","")
					if(status=="missing"): status += " &#9888;"

					tdEnthalpy = ""
					if(None in (productsEnthalpy+reactantsEnthalpy)):
						tdEnthalpy += "<td>missing enthalpy data"
					else:
						DeltaH = sum(productsEnthalpy)-sum(reactantsEnthalpy)
						DeltaH_K = kJmol2K*DeltaH
						tdEnthalpy += "<td>"+utils.htmlExp(DeltaH_K)
						if(DeltaH_K<0e0):
							tdEnthalpy += "<td>&#10004;&#10004;"
						elif(DeltaH_K>=0 and DeltaH_K<1e4):
							tdEnthalpy += "<td>&#10004;"
						else:
							tdEnthalpy += "<td>&#10006;"

					fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+row+"<td style=\"font-size:10px;\">"+status.upper()+tdEnthalpy+"\n")

			#dataBranch["presentBranches"] = []
			icount += 1
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

				#get DH for each reactant
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
				linkURL = "<a href=\"rate_"+myReaction.getReactionHash()+".html\">details</a> for "+myReaction.getVerbatimHtml()
				fout.write("<img src=\""+fnamePNG+"\" width=\"700px\" alt=\"&#9888; MISSING: "+myReaction.getVerbatim()+"\"><br>"+linkURL)
				fout.write("<tr height=\"10px\"><td>")
				icount += 1
			fout.write("</table>\n")

			#add footer
			fout.write(utils.getFooter("footer.php"))
			fout.close()

	#*******************
	#prepares HTML documentation for list of rate divided by number of intervals
	def makeHtmlMultipleRates(self,myOptions):

		#divide rates by number of intervals
		multipleRates = dict()
		for reaction in self.reactions:
			#get number of interval
			intervalsNumber = len(reaction.rate)
			#initialize the list for the corresponding dictionary key
			if(not(intervalsNumber in multipleRates)): multipleRates[intervalsNumber] = []
			#add the reaction to dict
			multipleRates[intervalsNumber].append(reaction)


		#prepare a file for each variable
		fname = "htmls/multipleRates.html"

		#standard header for rates
		tableHeader = "<tr>"+("<th>"*31)

		#open HTML file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<a href=\"index.html\">back</a><br>\n")
		#loop on multiple rates
		for intervalsNumber in sorted(multipleRates.keys())[::-1]:
			reactionBlock = multipleRates[intervalsNumber]
			fout.write("<br><br>\n")
			#write plural if necessary
			intervalString = "interval"+("s" if(intervalsNumber>1) else "")
			#table title
			fout.write("<p style=\"font-size:20px\">Reactions with "+str(intervalsNumber)+" "+intervalString+"</p>\n")
			#open reaction table
			fout.write("<table width=\"50%\">\n")
			fout.write(tableHeader+"\n")
			#reactions sorted by unsorted reaction hash (in the html list)
			reactionsSorted = sorted(reactionBlock, key=lambda x:x.getReactionHashUnsorted())

			icount = 0
			#loop on reactions
			for myReaction in reactionsSorted:
				bgcolor = ""
				if(icount%2!=0): bgcolor = utils.getHtmlProperty("tableRowBgcolor")
				fout.write("<tr valign=\"baseline\" bgcolor=\""+bgcolor+"\">"+myReaction.getReactionHtmlRow(mode="joints")+"\n")
				icount += 1
			fout.write(tableHeader+"\n")
			fout.write("</table>\n")


		#add footer
		fout.write(utils.getFooter("footer.php"))
		fout.close()


