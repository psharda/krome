import kexplorer_reaction,kexplorer_element
import sys,subprocess,os,glob,json,datetime

import itertools #added by Jels Boulangier 30/03/2017
########################################
#added by Jels Boulangier 05/04/2017
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

#plot arguments
font = {'size'   : 21}
lines = {'linewidth' : 5, 'markersize': 10, 'markeredgewidth': 3, }
savefig = {'dpi': 300, 'format': 'png', 'transparent': True,'bbox': 'tight'}
figure = {'figsize': (14, 10)}#,'autolayout':True}

plt.rc('font', **font)
plt.rc('lines', **lines)
plt.rc('savefig', **savefig)
plt.rc('figure', **figure)
plt.rc('xtick.major', size=10, width=1.5)
plt.rc('xtick.minor', size=5, width=1.5)
plt.rc('ytick.major', size=10, width=1.5)
plt.rc('ytick.minor', size=5, width=1.5)
########################################
class network:
	reactions = dict()
	elements = dict()#added by Jels Boulangier 05/04/2017

	minFlux = 1e-5 #minimum flux to plot
	plotLog = True #edge thickness is log of flux
	maxPenWidth = 5. #max edge thicknes, px
	xvarName = "" #name of independent variable
	xvarUnits = "" #units of independent variable
	xvarFormat = "%e" #format of independent variable
	#added by Jels Boulangier 05/04/2017
	minAbundance = 1e-10 #minimum mass fraction to plot
	maxAbundance = 1 #maximum mass fraction to plot


	#**********************
	#network contructor read kexplorer file
	def __init__(self,fileName,fileNameEvolution=None):

		#read data from file
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			#add reaction to this network object
			self.addReaction(srow)

		fh.close()

		#uses evolution data only if file name is present
		if(fileNameEvolution!=None):
			#added by Jels Boulangier 05/04/2017
			#read data from file
			fg = open(fileNameEvolution,"rb")
			#list of all elements
			fline = fg.readline().strip().split(" ")
			flineElem = fline[3:]
			#add element data per block in the inout file
			elemBlock = []
			for row in fg:
				srow = row.strip()
				if(srow==""):
					self.addElement(elemBlock,flineElem)
					elemBlock = []
				else:
					elemBlock.append(srow)

			fg.close()

			# determine (Tgas,xvar)-grid
			self.getRangeTgasXvar()




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

	#**********************
	#add element data per block in the input file
	#added by Jels Boulangier 05/04/2017
	def addElement(self,block,elemAll):
		newBlock = True
		#create element objects
		for el in elemAll:
			if(not(el in self.elements)):
				self.elements[el] = kexplorer_element.element()

		#add data for every element
		for row in block:
			arow = [x for x in row.split(" ") if(x!="")]
			(time, Tgas, xvar) = [float(x) for x in arow[0:3]]
			for el in elemAll:
				elIdx = elemAll.index(el)+3
				abundance = float(arow[elIdx])
				self.elements[el].addData(time, Tgas, xvar, abundance, newBlock)
			newBlock = False


	#******************
	#find all unique Tgas,xvar values
	#useful for determining the (tgas,xvar)-grid
	#added by Jels Boulangier 05/04/2017
	def getRangeTgasXvar(self):
		#find unique xvar/Tgas values ("H" is a random choice of element)
		xvarUniqueSet = set([i[0] for i in self.elements["H"].xvarData])
		tgasUniqueSet = set([i[0] for i in self.elements["H"].tgasData])
		#sorted list of unique xvar/Tgas values
		self.xvarUnique = sorted(list(xvarUniqueSet))
		self.tgasUnique = sorted(list(tgasUniqueSet))
		#range of xvar/Tgas values
		self.xvarRange = len(self.xvarUnique)
		self.tgasRange = len(self.tgasUnique)


	#******************
	#search for most fluxy reactions
	def findBest(self,threshold=.05,atom=None):

		#store atom base
		self.atomBase = atom

		self.bestReactions = []
		self.bestReactionsIdx = [] #added by Jels Boulangier 30/03/2017
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
				self.bestReactionsIdx.append(idx) #added by Jels Boulangier 30/03/2017

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
	#create network file with most fluxy reactions
	#added by Jels Boulangier 30/03/2017
	def networkBest(self,oldNetwork,newNetwork="BestNetwork.ntw"):

		def isa_group_separator(line):
		    return line=='\n'

		cnt = 0
		with open(oldNetwork, 'r') as fileInput, open(newNetwork, "w") as fileOutput:
			#split file into blocks separated by empty line, one block is one reaction
		    for key,group in itertools.groupby(fileInput,isa_group_separator):
		        if not key:
					#write most fluxy reactions in new network file (+ header)
		            if cnt in self.bestReactionsIdx or cnt==0:
		                reactionBlock = "".join(list(group))
		                fileOutput.write(reactionBlock + "\n")
		            cnt += 1

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
		for fname in glob.glob(pngFolder+"*.png"):
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
						if(not(product in graphMax[reactant])):
							graphMax[reactant][product] = 0e0
						#store maxium edge value
						graphMax[reactant][product] = \
							max(reaction.fluxNormTotData[i], graphMax[reactant][product])

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
						lineWidth = max(lineWidth-log10(self.minFlux),0e0) \
							* self.maxPenWidth/(-log10(self.minFlux))
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
					fout.write("\""+reactant+"\" -> \""+product+"\" "+style+";\n")
			fout.write("}\n")
			fout.close()

			#call graphviz to plot
			myCall = "dot tmp.dot -Tpng -o "+pngName
			subprocess.call(myCall.split(" "))

	#******************
	#make (T,xvar) color plot of ALL element abundances
	#added by Jels Boulangier 11/04/2017
	def abundanceColormapAll(self,elemInt=None,pngFolder="pngs"):
		#plot for all elements
		if not elemInt:
			for key in self.elements:
				self.abundanceColormap(key,pngFolder)
		#plot for elements of interest
		else:
			for elem in elemInt:
				self.abundanceColormap(elem,pngFolder)

	#******************
	#make (T,xvar) color plot of element abundances
	#added by Jels Boulangier 05/04/2017
	def abundanceColormap(self,atom,pngFolder="pngs"):

		print "Making abundance colormap of %s" %(atom)
		#get variable data after last time step
		x = [i[-1] for i in self.elements[atom].tgasData]
		y = [i[-1] for i in self.elements[atom].xvarData]
		z = [i[-1] for i in self.elements[atom].abundanceData]

		#create matrix for image plot
		Ncol = len(set(x))
		Nrow = len(set(y))
		z = np.reshape(z,(Nrow, Ncol))
		x = np.array(x)
		y = np.array(y)

		zMin = max(z.min(),self.minAbundance)
		zMax = min(z.max(),self.maxAbundance)
		zRange = zMax/zMin

		#do no make image if abundance is too small
		if zMax<self.minAbundance:
			print "Abundance of %s is zero everywhere" %(atom)
			return

		#create image
		plt.figure()
		if(zRange>10):
			#logaritmic colorbar
			plt.imshow(z, extent=(x.min(), x.max(),y.min(), y.max()), \
			interpolation='gaussian', cmap='afmhot',aspect='auto',origin='lower',\
			norm=colors.LogNorm(vmin=zMin, vmax=zMax))
		else:
			#linear colorbar
			plt.imshow(z, extent=(x.min(), x.max(),y.min(), y.max()), \
			interpolation='gaussian', cmap='afmhot',aspect='auto',origin='lower')
		#Acceptable interpolations are 'none', 'nearest', 'bilinear', 'bicubic',
		#'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser',
		#'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'

		#make plot labels
		plt.colorbar(label='Mass fraction')
		plt.yscale('log')
		plt.title('Fractional abundance of %s' %(atom))
		plt.xlabel('Temperature (K)')
		plt.ylabel(r'%s ($%s$)' %(self.xvarName,self.xvarUnits))
		#dump png file
		print "Dumping colormap of %s" %(atom)
		plt.savefig(pngFolder + '/%s' %(atom))
		plt.close()

	#******************
	#get the block index of a (Tgas,xvar) combo
	#added by Jels Boulangier 05/04/2017
	def getIdxTgasXvar(self,xvar,tgas):

		idxTgas = self.tgasUnique.index(tgas)
		idxXvar = self.xvarUnique.index(xvar)

		# TgasRange amount of blocks per xvar value
		idxBlock = idxXvar*self.tgasRange + idxTgas

		return idxBlock

	#******************
	#make plot of element evolution in (T,xvar)-space
	#added by Jels Boulangier 05/04/2017
	def abundanceEvolution(self,atom,tgasInt,xvarInt,pngFolder="pngs"):

		markers = ['D','o','*','+','s','x']
		colors = ['r','b','k','g','m','c','grey','brown','pink']
		markerHandles = []
		markerLabels = []
		colorHandles = []
		colorLabels = []

		#marker counter
		mIdx = 0
		#loop over all xvar
		for var in xvarInt:
			#marker
			m = markers[mIdx]
			#color counter
			cIdx = 0
			#loop over all temperatures
			for tgas in tgasInt:
				#color
				c = colors[cIdx]
				#find index of the wanted block of input file
				blockIdx = self.getIdxTgasXvar(var,tgas)

				#get abundance and time data
				x = self.elements[atom].timeData[blockIdx]
				y = self.elements[atom].abundanceData[blockIdx]
				plt.semilogy(x,y,'-'+m,color=c)

				#plot fake points for legend of colors
				if mIdx ==0:
					line, = plt.plot(-1,-1,'-',color=c)
					colorHandles.append(line)
					colorLabels.append("T = %i K" %(int(tgas)))

				cIdx += 1

			#plot fake points for legend of markers
			line, = plt.plot(-1,-1,'-'+m,color='grey')
			markerHandles.append(line)
			markerLabels.append(r"%s = %s $%s$" %(self.xvarName,str(var),self.xvarUnits))

			mIdx += 1

		#make legends for xvar and temperature
		legend1 = plt.legend(colorHandles,colorLabels,loc=2,\
		fancybox=True, shadow = True,bbox_to_anchor=(1, 1) )
		ax = plt.gca().add_artist(legend1)
		legend2 = plt.legend(markerHandles,markerLabels,loc=3,\
		fancybox=True, shadow=True, bbox_to_anchor=(1, 0))

		#make plot labels
		plt.title('Abundance evolution of %s' %(atom))
		plt.xlabel('Time (s)')
		plt.ylabel("Mass fraction")
		plt.xlim(xmin=0,xmax=x[-1])
		#dump png file
		print "Dumping abundance evolution of %s" %(atom)
		plt.savefig(pngFolder + '/ev%s' %(atom),bbox_extra_artists=[legend1,legend2])
		plt.close()


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
		hrefWiki = "<a href=\"https://bitbucket.org/tgrassi/krome/wiki/kexplorer\" " \
			+ " target=\"_blank\">kexplorer</a>"
		return "documentation generated with "+hrefWiki+" ("+bitbucketLanding+") - changeset: <a href=\""\
			+ bitbucketLink + "\" target=\"_blank\">" + changeset[:7] + "</a> - " + datenow

