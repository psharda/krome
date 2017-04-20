import kexplorer_reaction
import sys,subprocess,os,glob,json,datetime

class network:
	reactions = dict()

	minFlux = 1e-5 #minimum flux to plot
	plotLog = True #edge thickness is log of flux
	maxPenWidth = 5. #max edge thicknes, px
	xvarName = "" #name of independent variable
	xvarUnits = "" #units of independent variable
	xvarFormat = "%e" #format of independent variable

	#**********************
	#network contructor read kexplorer file
	def __init__(self,fileName):

		#read data from file
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			#add reaction to this network object
			self.addReaction(srow)

		fh.close()

	#**********************
	#add reaction to network using explore reaction file line format
	def addReaction(self,row):
		srow = row.strip()
		arow = [x for x in srow.split(" ") if(x!="")]
		idx = int(arow[0])
		#parse data
		(xvar, Tgas, flux, fluxNormMax, fluxNormTot) = [float(x) for x in arow[1:6]]

		#parse verbatim reaction and create reaction object
		if(not(idx in self.reactions)):
			rname = (" ".join(arow[6:]))
			self.reactions[idx] = kexplorer_reaction.reaction(rname)

		#add data to reaction
		self.reactions[idx].addData(xvar, Tgas, flux, fluxNormMax, fluxNormTot)


	#******************
	#search for most fluxy reactions
	def findBest(self,threshold=.05,atom=None):

		#store atom base
		self.atomBase = atom

		self.bestReactions = []
		#loop on reactions to find criteria
		for (idx,reaction) in self.reactions.iteritems():
			if(atom!=None):
				hasAtom = False
				for species in (reaction.reactants+reaction.products):
					if(atom in species): hasAtom = True
				if(not(hasAtom)): continue
			#if percent flux greater than treshold then store reaction
			if(max(reaction.fluxNormTotData)>threshold):
				self.bestReactions.append(reaction)

		#if no reactions found rise error
		if(len(self.bestReactions)==0):
			print "ERROR: no reactions found with threshold "+str(threshold)
			print " and atom criterium ",atom
			sys.exit()
		print "reactions found ",len(self.bestReactions)

	#*****************
	#list most fluxy reactions on screen
	def listBest(self):
		for reaction in self.bestReactions:
			print reaction.verbatim

	#****************
	#buld a graph map with the list of reaction partners
	def buildGraphMap(self,reactions):
		graphMap = dict()
		#loop on reactions
		for reaction in reactions:
			#loop on reactants
			for ridx in range(len(reaction.reactants)):
				reactant = reaction.reactants[ridx]
				#init graph dict
				if(not(reactant in graphMap)):
					graphMap[reactant] = dict()
				#loop on reactions
				for product in reaction.products:
					#init partners list
					if(not(product in graphMap[reactant])):
						graphMap[reactant][product] = []
					#add reaction partners
					partners = [reaction.reactants[i] for i in range(len(reaction.reactants)) if(i!=ridx)]
					#add partners
					graphMap[reactant][product] += partners

		#graphMap as class attribute
		self.graphMap = graphMap

	#****************
	#dump data to JSON
	def dumpJSON(self):

		jsonData = []
		#get dependent variable data
		xvarData = self.bestReactions[0].xvarData
		#count data
		numberData = len(xvarData)
		for i in range(numberData):
			stepData = []
			for reaction in self.bestReactions:
				rData = {"verbatim":reaction.verbatim, \
					"htmlVerbatim":reaction.htmlVerbatim, \
					"flux": reaction.fluxData[i], \
					"fluxNormTot":reaction.fluxNormTotData[i]}
				stepData.append(rData)
			stepData = sorted(stepData,key=lambda x:x["flux"])[::-1]
			jsonData.append({"data":stepData,"xvar":xvarData[i]})

		with open('data.json', 'w') as outfile:
			json.dump(jsonData, outfile)

	#****************
	#dump graph to png
	def dumpBest(self,pngFolder="pngs/"):
		from math import log10


		#if png folder not present than create
		if(not(os.path.exists(pngFolder))):
			os.makedirs(pngFolder)

		#clear pngs
		for fname in glob.glob("pngs/*.png"):
			os.remove(fname)

		#get partners
		self.buildGraphMap(self.bestReactions)

		#get dependent variable data
		xvarData = self.bestReactions[0].xvarData
		#count data
		numberData = len(xvarData)

		self.dumpJSON()
		self.makeHTML(numberData)

		#loop on data
		for i in range(numberData):
			#png filename
			pngName = pngFolder+"/graph_"+str(int(1e6+i))+".png"
			pngName = pngName.replace("//","/")
			print "generating "+str(i+1)+"/"+str(numberData)+" as "+pngName
			graphMax = dict()
			allSpecies = []
			#loop on reactions
			for reaction in self.bestReactions:
				#loop on reactants
				for reactant in reaction.reactants:
					#store all species
					allSpecies.append(reactant)
					#init dictionary of products
					if(not(reactant in graphMax)): graphMax[reactant] = dict()
					#loop on products
					for product in reaction.products:
						#store all species
						allSpecies.append(product)
						#init variable data
						if(not(product in graphMax[reactant])): graphMax[reactant][product] = 0e0
						#store maxium edge value
						graphMax[reactant][product] = max(reaction.fluxNormTotData[i],graphMax[reactant][product])

			#create temporary dot file
			fout = open("tmp.dot","w")
			fout.write("digraph G{\n")
			xvarName = ""
			if(self.xvarName!=""): xvarName = self.xvarName+" = "
			xvar = xvarName+(self.xvarFormat % xvarData[i])+" "+self.xvarUnits
			tgasVar = "Tgas = "+str(round(reaction.tgasData[i],1))+" K"
			#add independent variable to PNG
			fout.write("label=\""+xvar+", "+tgasVar+"\";\n")
			#unique list of species
			allSpecies = list(set(allSpecies))
			#loop on species to have circle nodes
			for species in allSpecies:
				#skip species if atomBase is present
				if(self.atomBase!=None):
					if(not(self.atomBase in species)): continue
				fout.write("\""+species+"\" [shape=circle];\n")

			#loop on graph map, reactants
			for (reactant,products) in graphMax.iteritems():
				#loop on graph map, products
				for (product,maxval) in products.iteritems():

					#check if reactant+prducts has atom base, if not skip
					if(self.atomBase!=None):
						reactantHasBase = (self.atomBase in reactant)
						productHasBase = (self.atomBase in product)
						if(not(reactantHasBase) or not(productHasBase)): continue

					if(self.plotLog):
						#get log of flux
						lineWidth = log10(graphMax[reactant][product]+1e-40)
						#set a limit to flux size
						lineWidth = max(lineWidth-log10(self.minFlux),0e0)*self.maxPenWidth/(-log10(self.minFlux))
					else:
						lineWidth = graphMax[reactant][product]
						if(lineWidth<self.minFlux): lineWidth = 0e0
						lineWidth *= self.maxPenWidth

					#default color is black
					color = "#000000"
					#when below threshold draw light blue edge
					if(lineWidth==0e0):
						lineWidth = 1
						color = "#9ACEEB"
					#edge style
					style = "[penwidth="+str(lineWidth)+", color=\""+color+"\", arrowsize=.5]"
					#write to dot file
					fout.write("\""+reactant+"\" -> \""+product+"\" "+style+";\n") #self.graphMap[reactant][product]
			fout.write("}\n")
			fout.close()

			#call graphviz to plot
			myCall = "dot tmp.dot -Tpng -o "+pngName
			subprocess.call(myCall.split(" "))

	#**************
	def loadHeader(self,fname,maxSteps):
		fh = open(fname,"rb")
		allFile = ""
		for row in fh:
			allFile += row.replace("#MAXFILE#",str(maxSteps))
		fh.close()
		return allFile

	#**************
	def loadFooter(self,fname):
		fh = open(fname,"rb")
		allFile = ""
		for row in fh:
			allFile += row.replace("#FOOTER_INFO#",self.getFooterInfo())
		fh.close()
		return allFile

	#**************
	def makeHTML(self,maxSteps):
		fh = open("kexplorer.php","rb")
		fout = open("index.html","w")
		for row in fh:
			row = row.replace("#HEADER#",self.loadHeader("header.php",maxSteps))
			row = row.replace("#FOOTER#",self.loadFooter("footer.php"))
			fout.write(row)
		fh.close()
		fout.close()

	#***********************
	def getFooterInfo(self):
		#name of the git master file
		masterfile = "../../.git/refs/heads/master"
		changeset = ("x"*7) #default unknown changeset
		#if git master file exists grep the changeset
		if(os.path.isfile(masterfile)):
			changeset = open(masterfile,"rb").read()

		datenow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		bitbucketLink = "https://bitbucket.org/tgrassi/krome/commits/"+changeset
		bitbucketLanding = "<a href=\"https://bitbucket.org/tgrassi/krome\" target=\"_blank\">KROME</a>"
		hrefWiki = "<a href=\"https://bitbucket.org/tgrassi/krome/wiki/kexplorer\" target=\"_blank\">kexplorer</a>"
		return "documentation generated with "+hrefWiki+" ("+bitbucketLanding+") - changeset: <a href=\""\
			+ bitbucketLink + "\" target=\"_blank\">" + changeset[:7] + "</a> - " + datenow



