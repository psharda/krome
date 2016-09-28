import sys,glob,os
import reaction,utils,options

class network:

	species = None

	#************************
	#class constructor
	def __init__(self,myOptions):

		#default format
		reactionFormat = "@idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"

		#get network file name
		fileName = myOptions.network

		#read atoms
		atomSet = utils.getAtomSet("atomlist.dat")
		#read shortcuts
		shortcuts = utils.getShortcuts()

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

			#skip other tokens
			if(srow.startswith("@")): continue

			#parse row line for reaction
			myReaction = reaction.reaction(srow,reactionFormat,atomSet)

			#add parsed reaction to reactions structure in network
			self.reactions.append(myReaction)
		fh.close()

		#merge reactions with multiple rates
		print "merging multiple reactions"
		self.mergeReactions()

		#get ranges from option file
		varRanges = dict()
		for rng in myOptions.range:
			#store range name and range limits
			(rangeName,rangeValue) = [x.strip() for x in rng.split("=")]
			varRanges[rangeName] = [float(x) for x in rangeValue.split(",")]
			print rangeName, varRanges[rangeName]

		#clear folders (html and png)
		self.clearFolders()

		#create reaction index page
		self.makeHtmlReactionIndex()

		#prepare html pages for species
		for mySpecies in self.getSpecies():
			mySpecies.makeHtmlPage(self)

		#plotting rates
		print "plotting rates..."
		icount = 0
		#loop on reactions to evalute and plot
		for myReaction in self.reactions:
			print str(int(icount*1e2/len(self.reactions)))+"%", myReaction.getVerbatim()
			myReaction.plotRate(shortcuts,varRanges)
			myReaction.makeHtmlPage(myOptions)
			icount += 1

	#**************
	def getSpecies(self):
		if(self.species!=None): return self.species
		self.species = []
		for reaction in self.reactions:
			self.species += reaction.reactants + reaction.products

		self.species = list(set(self.species))
		return self.species

	#**************
	#clean temporary folders
	def clearFolders(self):
		#files to be deleted
		folders = {"pngs":"png", \
			"htmls":"html"}

		#loop on folders and extensions
		for (path,extension) in folders.iteritems():
			print "removing "+extension+"s..."
			#get files list
			filelist = glob.glob(path+"/*."+extension)
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


	#****************
	#create reaction list index as html page
	def makeHtmlReactionIndex(self):

		fname = "htmls/indexReactions.html"

		#open file to write
		fout = open(fname,"w")
		#add header
		fout.write(utils.getFile("header.php"))
		fout.write("<br><br>\n")
		#reaction table
		fout.write("<table width=\"50%\">\n")
		fout.write("<tr><th><th><th>\n")
		#loop on reactions
		for myReaction in self.reactions:
			fout.write("<tr>"+myReaction.getReactionHtmlRow()+"\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("</table>\n")

		#add footer
		fout.write(utils.getFile("footer.php"))
		fout.close()


