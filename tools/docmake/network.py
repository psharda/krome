import sys,glob,os
import reaction,utils,options

class network:

	#************************
	#class constructor
	def __init__(self,myOptions):

		#default format
		reactionFormat = "@idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"

		fileName = myOptions.network

		atomSet = utils.getAtomSet("atomlist.dat")
		shortcuts = utils.getShortcuts()

		self.reactions = []
		print "reading network "+fileName
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue

			if(srow.startswith("@format:")):
				reactionFormat = srow
				continue

			if(srow.startswith("@var:")):
				(variable,expression) = [x.strip() for x in srow.replace("@var:","").split("=")]
				shortcuts[variable] = expression

			if(srow.startswith("@")): continue
			myReaction = reaction.reaction(srow,reactionFormat,atomSet)
			self.reactions.append(myReaction)
		fh.close()

		print "merging multiple reactions"
		self.mergeReactions()

		#get ranges from option file
		varRanges = dict()
		for rng in myOptions.range:
			(rangeName,rangeValue) = [x.strip() for x in rng.split("=")]
			varRanges[rangeName] = [float(x) for x in rangeValue.split(",")]
			print rangeName, varRanges[rangeName]

		#clear folders
		self.clearFolders()

		#create reaction index page
		self.makeHtmlReactionIndex()

		#plotting rates
		print "plotting rates..."
		for myReaction in self.reactions:
			print myReaction.getVerbatim()
			myReaction.plotRate(shortcuts,varRanges)
			myReaction.makeHtmlPage(myOptions)

	#**************
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
	def makeHtmlReactionIndex(self):

		fname = "htmls/indexReactions.html"

		fout = open(fname,"w")
		fout.write(utils.getFile("header.php"))
		fout.write("<br><br>\n")
		fout.write("<table width=\"50%\">\n")
		fout.write("<tr><th><th><th>\n")
		for myReaction in self.reactions:
			fout.write("<tr>"+myReaction.getReactionHtmlRow()+"\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("</table>\n")

		fout.write(utils.getFile("footer.php"))
		fout.close()


