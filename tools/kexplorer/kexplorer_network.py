import kexplorer_reaction,kexplorer_element,kexplorer_utils, figureSettings
import sys,subprocess,os,glob,json,datetime

import itertools #added by Jels Boulangier 30/03/2017
########################################
#added by Jels Boulangier 05/04/2017
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

########################################
class network:
	reactions = dict()

	minFlux = 1e-5 #minimum flux to plot
	plotLog = True #edge thickness is log of flux
	maxPenWidth = 5. #max edge thicknes, px
	xvarName = "" #name of independent variable
	xvarUnits = "" #units of independent variable
	xvarFormat = "%e" #format of independent variable
	#added by Jels Boulangier 05/04/2017
	minAbundance = 1e-20 #minimum mass fraction to plot
	maxAbundance = 1e10 #maximum mass fraction to plot
	reaFormat = "idx,R,R,R,P,P,P,P,Tmin,Tmax,rate" #default format
	rateLength = 100 #maximum length of rate

	#**********************
	#network contructor read kexplorer file
	def __init__(self,fileName,fileNameEvolution=None,elemInterest=None):

		#make dict to store element objects
		self.elements = dict()
		#read data from file
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			#add reaction to this network object
			self.addReaction(srow)

		fh.close()

		#uses evolution data only if file name is present
		if fileNameEvolution:
			#added by Jels Boulangier 05/04/2017
			#read data from file
			fg = open(fileNameEvolution,"rb")
			#list of all elements
			self.allSpecieslist = fg.readline().strip().split(" ")[3:]
			if elemInterest:
				flineElem = elemInterest
			else:
				flineElem = self.allSpecieslist
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
				elIdx = self.allSpecieslist.index(el)+3
				abundance = float(arow[elIdx])
				self.elements[el].addData(time, Tgas, xvar, abundance, newBlock)
			newBlock = False


	#******************
	#find all unique Tgas,xvar values
	#useful for determining the (tgas,xvar)-grid
	#added by Jels Boulangier 05/04/2017
	def getRangeTgasXvar(self):
		#find unique xvar/Tgas values (take first as random choice of element)
		xvarUniqueSet = set([i[0] for i in self.elements[self.elements.keys()[0]].xvarData])
		tgasUniqueSet = set([i[0] for i in self.elements[self.elements.keys()[0]].tgasData])
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
	def abundanceColormapAll(self,elemInt=None,timeEvolution=False,pngFolder="pngs"):

		if elemInt:
			#do for elements of interest
			speciesTodo = elemInt
		else:
			#do for all elements
			speciesTodo = self.elements

		if timeEvolution:
			# number of time grid points
			timeStart = 0
			timeStop = len(self.elements[self.elements.keys()[0]].timeData[0])
		else:
			timeStart = -1
			timeStop = 0
		for species in speciesTodo:
			for timeIndex in range(timeStart,timeStop):
				if timeIndex==timeStart:
					globalMinimum = 1e64
					globalMaxmimum = 0
					for i in self.elements[species].abundanceData:
						for j in i:
							if j < globalMinimum and j != 0:
								globalMinimum = j
							if j > globalMaxmimum:
								globalMaxmimum = j

					limits = [globalMinimum, globalMaxmimum]

				self.abundanceColormap(species, timeIndex, limits, pngFolder)


	#******************
	# Make a video of all the abundance colormap to visualise the
	# temporal evolution in temperature-xvar space
	# It reads in png files that were created with
	# abundanceColormapAll(timeEvolution=True)
	def makeEvolutionVideo(self, elemInt=None,pngFolder="pngs"):

		if elemInt:
			#do for elements of interest
			speciesTodo = elemInt
		else:
			#do for all elements
			speciesTodo = self.elements

		for species in speciesTodo:
			# change video setting if desired
			# this uses ffmpeg
			os.system("ffmpeg -framerate 3 -pattern_type glob -i "
					"\'" + pngFolder + "/" + species + "_*.png\' "
					"-c:v libx264 -s 1920x1080 -r 30 -pix_fmt yuv420p "
					+ pngFolder + "/" + species + "evolution.mp4" )


	#******************
	#make (T,xvar) color plot of element abundances
	#added by Jels Boulangier 05/04/2017
	def abundanceColormap(self,atom,idxTime=-1,limits=None,pngFolder="pngs"):

		print "Making abundance colormap of %s" %(atom)

		if idxTime!=-1:
			evolution = True
		else:
			evolution = False


		#get variable data after last time step
		x = [i[idxTime] for i in self.elements[atom].tgasData]
		y = [i[idxTime] for i in self.elements[atom].xvarData]
		z = [i[idxTime] for i in self.elements[atom].abundanceData]

		#create matrix for image plot
		Ncol = len(set(x))
		Nrow = len(set(y))
		z = np.reshape(z,(Nrow, Ncol))
		x = np.reshape(x,(Nrow, Ncol))
		y = np.reshape(y,(Nrow, Ncol))

		if evolution:
			zMin = max(limits[0],self.minAbundance)
			zMax = min(limits[1],self.maxAbundance)

		else:
			zMin = max(z.min(),self.minAbundance)
			zMax = min(z.max(),self.maxAbundance)
		zRange = zMax/zMin

		#do no make image if abundance is too small
		if zMax<self.minAbundance:
			print "Abundance of %s is zero everywhere" %(atom)
			return

		# replace all extremely small numbers with a 'okay' small number
		# to avoid NaNs
		z[z < 1e-60] = 1e-60

		#create image
		if(zRange>10):
			# logaritmic colorbar
			# Two options are presented here, and the user can choose
			# by commenting one or the other

			### Color plot with a gradient color bar
			# rasterized=True is needed for pdfs.

			# plt.pcolormesh(x, y, z, cmap='viridis', rasterized=True,
			# norm=colors.LogNorm(vmin=zMin, vmax=zMax))

			### Color plot with filled contours and (connected) contour lines
			# NOTE: that there is a bug in matplotlib when using 'extend':
			# "ValueError: extend kwarg does not work yet with log scale"
			# Note that there is no 'extend' arrow on the color bar...
			#
			## OPTION 1: do a quick and dirty fix by replacing all
			# values below lower limit, with lower limit.
			#
			## OPTION 2: wait till my pull request 10705 gets accepted
			# (https://github.com/matplotlib/matplotlib/pull/10705)
			# OR fix this in your own matplotlib version.
			# e.g. path where the file typically resides:
			# /anaconda3/envs/py27/lib/python2.7/site-packages/matplotlib/contour.py
			# Step 1: in def _init_, remove the ValueError if statement
			# Step 2: in def _process_levels,
			# change :
			# if self.extend in ('both', 'min'):
			# 	self.layers[0] = -1e150
			# to:
			# if self.extend in ('both', 'min'):
			# 	if self.logscale:
            #     self.layers[0] = 1e-150
            # 	else:
            #     self.layers[0] = -1e150
			# Reason: very large negative number gives NaN when using log,
			# so replace with very small number.
			#
			## OPTION 1:

			# exp_min = np.floor(np.log10(zMin)-1)
			# exp_max = np.ceil(np.log10(zMax)+1)
			# lev_exp = np.arange(exp_min, exp_max)
			# levs = np.power(10, lev_exp)
			# #replace all values below lower limit, with lower limit.
			# z[z < levs[0] ] = levs[0]
			#
			# plt.contourf(x, y, z, levels=levs,
			# 			cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

			## OPTION 2:
			exp_min = np.floor(np.log10(zMin))
			exp_max = np.ceil(np.log10(zMax)+1)
			lev_exp = np.arange(exp_min, exp_max)
			levs = np.power(10, lev_exp)
			# Try best/newest matplotlib option (OPTION 2).
	        # If not present, then do OPTION 1.
	        try:
	            plt.contourf(x, y, z, levels=levs, extend='both',
	                        cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

	        except(ValueError,), err:
	            # replace this with the pcolormesh option if desired.
	            exp_min = np.floor(np.log10(zMin)-1)
	            exp_max = np.ceil(np.log10(zMax)+1)
	            lev_exp = np.arange(exp_min, exp_max)
	            levs = np.power(10, lev_exp)
	            #replace all values below lower limit, with lower limit.
	            z[z < levs[0] ] = levs[0]

	            plt.contourf(x, y, z, levels=levs,
	            			cmap='viridis', norm=colors.LogNorm(vmin=zMin, vmax=zMax))

		else:
			#linear colorbar
			# same as in the logaritmic case:

			# plt.pcolormesh(x, y, z, cmap='viridis', rasterized=True)

			plt.contourf(x, y, z, cmap='viridis')

		#make plot labels
		#BUG sometime there are no labels/ticks on the colorbar...
		plt.colorbar(label='Number density', extend='min')
		plt.yscale('log')
		if evolution:
			plt.title('Abundance of %s time step %03i' %(atom,idxTime))
		else:
			plt.title('Fractional abundance of %s' %(atom))
		plt.xlabel('Temperature (K)')
		plt.ylabel(r'%s (%s)' %(self.xvarName,self.xvarUnits))
		#dump png file
		print "Dumping colormap of %s" %(atom)
		if evolution:
			plt.savefig(pngFolder + '/%s_%03i.png' %(atom,idxTime),
						format='png', dpi=200)
		else:
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
	def abundanceEvolution(self,species,tgasInt,xvarInt,pngFolder="pngs"):

		markers = ['D','+','*','o','s','x','v','3','d','8']
		colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
		styles = ["-",":","--","-."]
		markerHandles = []
		markerLabels = []
		styleHandles = []
		styleLabels = []
		colorHandles = []
		colorLabels = []

		#style counter
		sIdx = 0
		#loop over all xvar
		for var in xvarInt:
			#style
			s = styles[sIdx]
			#color counter
			mIdx = 0
			#loop over all temperatures
			for tgas in tgasInt:
				#color
				m = markers[mIdx]
				#find index of the wanted block of input file
				blockIdx = self.getIdxTgasXvar(var,tgas)
				#marker counter
				cIdx = 0
				#loop over al species
				for spec in species:
					#marker
					c = colors[cIdx]
					#get abundance and time data
					x = self.elements[spec].timeData[blockIdx]
					y = self.elements[spec].abundanceData[blockIdx]
					plt.semilogy(x,y,s+m,color=c,markevery=5)

					#plot fake points for legend of markers
					if sIdx ==0 and mIdx ==0:
						line, = plt.plot(-1,-1,'-',color=c)
						colorHandles.append(line)
						colorLabels.append("%s" %(spec))

					cIdx += 1

				#plot fake points for legend of colors
				if sIdx ==0:
					line, = plt.plot(-1,-1,m,color='grey')
					markerHandles.append(line)
					markerLabels.append("T = %i K" %(int(tgas)))

				mIdx += 1

			#plot fake points for legend of styles
			line, = plt.plot(-1,-1,s,color='grey')
			styleHandles.append(line)
			styleLabels.append(r"%s = %s %s" %(self.xvarName,str(var),self.xvarUnits))

			sIdx += 1

		#make legends for xvar and temperature
		legend1 = plt.legend(colorHandles,colorLabels,loc=2,
		fancybox=True, shadow = True,bbox_to_anchor=(1, 1) )
		ax = plt.gca().add_artist(legend1)
		legend2 = plt.legend(markerHandles,markerLabels,loc=2,
		fancybox=True, shadow=True, bbox_to_anchor=(1, 0.7))
		ax = plt.gca().add_artist(legend2)
		legend3 = plt.legend(styleHandles,styleLabels,loc=3,
		fancybox=True, shadow=True, bbox_to_anchor=(1, 0))

		#make plot labels
		plt.title('Abundance evolution of %s' %(spec))
		plt.xlabel('Time (s)')
		plt.ylabel("Mass fraction")
		plt.xlim(xmin=0,xmax=x[-1])
		plt.ylim(ymin=self.minAbundance)
		#dump png file
		print "Dumping abundance evolution of %s" %(spec)
		# plt.show()
		plt.savefig(pngFolder + '/ev%s' %(spec),bbox_extra_artists=[legend1,legend2,legend3])
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
