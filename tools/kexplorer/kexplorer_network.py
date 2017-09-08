import kexplorer_reaction,kexplorer_element,kexplorer_utils
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
	reaFormat = "idx,R,R,R,P,P,P,P,Tmin,Tmax,rate" #default format

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

		cnt = 0
		with open(oldNetwork, 'r') as fileInput, open(newNetwork, "w") as fileOutput:
			#split file into blocks separated by empty line, one block is one reaction
		    for key,group in itertools.groupby(fileInput,kexplorer_utils.blockSeparator):
		        if not key:
					#write most fluxy reactions in new network file (+ header)
		            if cnt in self.bestReactionsIdx or cnt==0:
		                reactionBlock = "".join(list(group))
		                fileOutput.write(reactionBlock + "\n")
		            cnt += 1

	#****************
	#create latex format for network
	def network2latex(self,networkFile,networkLatex="NetworkLatex.log" ):

		reactionFormat = self.reaFormat #use default reaction format
		reactionCnt = 0 #reaction counter

		with open(networkFile, 'r') as fileInput, open(networkLatex, "w") as fileOutput:
			#dump header of the file
			self.dumpLatexTableHeader(fileOutput)
			#split file into blocks separated by empty line, one block is one reaction
			for key,group in itertools.groupby(fileInput,kexplorer_utils.blockSeparator):
				if not key:
					newReaction = True #keep track if new reaction or not (double rates)
					reactionBlockList = list(group)
					varList = [] #dictionary with possibly needed rate variables
					for item in reactionBlockList:
						if item.startswith("#"): continue
						if item.startswith("@var"):
							varLine  = item.split(":")[-1].split("=")
							varValue = varLine[-1].strip() #variable value
							varName  = varLine[0].strip()  #variable

							#adapt string to avoid math mistakes later on
							varValue = "(" + varValue + ")"
							varList.append((varName,varValue))

						#use new format if detected, if not, keep using previous
						if item.startswith("@format"):
							reactionFormat = item.split(":")[-1]

						#skip other special lines
						elif item.startswith("@"):
							continue

						elif int(item[0]):
							reactionInfo = item.split(",")
							#get reaction rare
							reactionRate = reactionInfo[-1].strip("\n")
							#transform into LaTeX format
							reactionRateTex = self.rate2latex(reactionRate,varList)

							if newReaction:
								formatList = reactionFormat.split(",")
								firstReactantIdx = kexplorer_utils.indicesElemList(formatList,"R")[0]
								lastProductIdx = kexplorer_utils.indicesElemList(formatList,"P")[-1]

								#increment reaction counter
								reactionCnt += 1
								reactionIdx = str(reactionCnt)

								formatReaProd = formatList[firstReactantIdx:lastProductIdx+1]
								verbReaction = reactionInfo[firstReactantIdx:lastProductIdx+1]
								verbReactionTex = self.reaction2latex(verbReaction,formatReaProd)

							else:
								verbReactionTex = ""
								reactionIdx = ""
							#dump line to output file
							#fileOutput.write(reactionRateTex + "\n")
							#fileOutput.write(verbReactionTex + "\n")
							newReaction = False
							#temporatry empty
							tempRangeTex = ""
							refTex = ""

							latexColums = [reactionIdx,verbReactionTex,reactionRateTex,tempRangeTex,refTex]
							self.dumpLatexTable(latexColums,fileOutput)

	#****************
	#KROME reaction format to LaTeX format
	def reaction2latex(self,reaction,format):

		reaTex = ""
		totalReact = format.count("R")
		totalProd  = format.count("P")

		#***************-
		def reaProdTex(out,it):
			elem = reaction[it]
			#skip empty elem (due to non-asdapted format)
			if elem =="": return ""
			#sub and superscrits to LaTeX format
			elem = kexplorer_utils.subSuper2latex(elem)
			elem = kexplorer_utils.special2latex(elem)
			if it==0 or it==totalReact: out = elem
			else: out = " $+$ " + elem
			return out

		#loop over all reactants
		for i in range(totalReact): reaTex += reaProdTex(reaTex,i)
		if totalReact==1: reaTex += " $+$ CR"
		# add reaction arrow
		reaTex = reaTex + " $\\to$ "
		#loop over all products
		for i in range(totalReact,totalReact+totalProd): reaTex += reaProdTex(reaTex,i)
		#add photon to reaction with only one product
		if totalProd==1: reaTex += " $+$ $\\gamma$ "

		return reaTex

	#****************
	#KROME network format to LaTeX format
	def rate2latex(self,rate,varList):
		import re
		import sympy as sp

		#list of symbols you want to keep in the LaTeX format
		Tsymbols = ["T","(T/300)", "T_{e}"]
		T, T32, Te = sp.symbols(Tsymbols)
		exp = sp.Symbol("exp")
		ln = sp.Symbol("ln")
		log = sp.Symbol("log")
		sqrt = sp.Symbol("sqrt")

		#put all variables with corresponding values in rate
		if varList:
			#loop needs to be reversed order for variable dependencies
			for var in reversed(varList):
				if var[0] in rate:
					rate = rate.replace(var[0],var[1])

		#list with all temperature shortcuts element = (var, replaceWith)
		tempShort = kexplorer_utils.getShortcuts()
		#replace shortcuts, loop needs to be reversed order for variable dependencies
		#skip T32, to keep as symbol
		for tshort in reversed(tempShort[2:]):
			rate = rate.replace(tshort[0],tshort[1])

		#make sympy friendly
		rate = rate.replace("d","e")
		rate = rate.replace("log", "ln")
		rate = rate.replace("ln10", "log")

		#transform to LaTeX format
		#keep trying i
		while True:
			try:
				rateTex = sp.latex(eval(rate))
				break
			#undefined variable will become a symbol
			except (NameError,),err:
				print "Name error in rate", err
				varIssue = str(err).split("'")[1]
				symb = varIssue + " = sp.Symbol(\""+varIssue+"\")"
				exec(symb)

			#special case rate	will be prited as it is
			except (SyntaxError,),err:
				print "Syntax Error in rate", err
				return rate

		#fix mistakes by sympy
		#no 10^{} for short rates
		for t in Tsymbols:
			if t in rateTex:
				break
			else:
				rateTex = rateTex.replace("e-","\\cdot 10^{-") + "}"
				break

		#it automatically makes a fraction out of negative exponents
		stringFrac = r'\frac{1}{'
		if stringFrac in rateTex:
			for Tsym in Tsymbols:
				rateTex = rateTex.replace(stringFrac + Tsym + "^{", "{"+Tsym+"^{-")
			#change the order of the factors to match modified Arrhenius
			#avoid chaning stuff in more complex rates
			if len(rateTex) < 100:
				rateTexSplit = re.split("\}\}",rateTex)
				beta = rateTexSplit[0][1:]+"}" #beta containing factor
				if "exp" in rateTexSplit[-1]:
					alphaGamma = rateTexSplit[-1].split("\operatorname") #alpha and gamma factor
					rateTex = alphaGamma[0] + beta + "\operatorname" + alphaGamma[1]
				else:
					rateTex = rateTexSplit[-1] + beta

		#mistake: large fractions
		#solution: to the power -1 (solution can be improved)
		if rateTex.startswith(r"\frac"):
			cnt = 0
			pieces = []
			#get content between parenteses for each level
			parenticList = kexplorer_utils.getParentheticContents(rateTex)
			#only get the first \frac{}{} parts
			for part in parenticList:
				if part[0]==0 and cnt<3:
					pieces.append(part[1])
					cnt = cnt +1

			#replace \frac{a}{b} with a*b^{-1}
			firstTerm = pieces[0]
			secondTerm = pieces[1]

			fracStringOriginal = r"\frac{"+firstTerm+"}{"+secondTerm+"}"

			# Put parenteses around first term if composite term
			if " + " not in firstTerm or " - " not in firstTerm:
				fracStringReplace = firstTerm+"\left["+secondTerm+r"\right]^{-1}"
			else:
				fracStringReplace = "\left["+firstTerm+r"\right]\left["+secondTerm+r"\right]^{-1}"

			rateTex = rateTex.replace(fracStringOriginal,fracStringReplace)
		#add LaTeX symbols
		rateTex = "$" + rateTex + "$"
		return rateTex

	#****************
	#dump colums in LateX table format
	def dumpLatexTableHeader(self,tableFile="latexTable.log"):

		tableFile.write("#********************************\n")
		tableFile.write("#Chemical reaction network (LaTeX format)\n")
		tableFile.write("#Table columns format {"+5*"l"+"}\n\n")

	#****************
	#build colums in LateX table format
	def dumpLatexTable(self,colums,tableFile="latexTable.log"):
		row = ""
		cnt = 1
		for el in colums:
			if cnt==len(colums): row += el + " \\\\"
			else: row += el + " & "
			cnt += 1

		tableFile.write(row + "\n")

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
