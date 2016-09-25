import sys,species,utils,os
from math import log10,log,exp,sqrt
class reaction:

	index = -1
	verbatim = None
	reactants = []
	products = []
	Tmin = []
	Tmax = []
	rate = []
	verbatim = verbatimLatex = reactionHash = None
	nameHtml = verbatimHtml = reactionHtmlRow = None
	evaluation = []

	#********************
	#parse csv reaction
	def __init__(self,row,reactionFormat,atomSet):
		if(not(reactionFormat.startswith("@format:"))):
			sys.exit("ERROR: wrong format "+reactionFromat)
		splitFormat = reactionFormat.replace("@format:","").split(",")
		splitFormat = [x.lower() for x in splitFormat]
		arow = row.strip().split(",",len(splitFormat))
		arow = [x.strip() for x in arow]

		specials = ["","G","CR"]

		self.rate = []
		self.reactants = []
		self.products = []
		self.Tmin = [None]
		self.Tmax = [None]

		for i in range(len(splitFormat)):
			part = splitFormat[i]
			if(part=="idx"):
				self.index = int(arow[i])
			elif(part=="r"):
				if(arow[i].upper() in specials): continue
				spec = species.species(arow[i],atomSet)
				self.reactants.append(spec)
			elif(part=="p"):
				if(arow[i].upper() in specials): continue
				spec = species.species(arow[i],atomSet)
				self.products.append(spec)
			elif(part=="tmin"):
				self.Tmin = [arow[i].replace("d","e")]
				if(arow[i].lower()=="none" or arow[i]==""): self.Tmin = [None]
			elif(part=="tmax"):
				self.Tmax = [arow[i].replace("d","e")]
				if(arow[i].lower()=="none" or arow[i]==""): self.Tmax = [None]
			elif(part=="rate"):
				self.rate.append(arow[i])
			else:
				print "ERROR: unknow format element "+part
				sys.exit(reactionFormat)
		self.check()

	#**************
	def check(self):
		reactantsCharge = sum([x.charge for x in self.reactants])
		productsCharge = sum([x.charge for x in self.products])
		if(reactantsCharge!=productsCharge):
			print "ERROR: charge problems"
			print self.getVerbatim()
			sys.exit()

		reactantsMass = sum([x.mass for x in self.reactants])
		productsMass = sum([x.mass for x in self.products])
		me = 9.10938356e-28 #electron mass, g
		if(abs(reactantsMass-productsMass)>=me/2.):
			print "ERROR: mass problems"
			print self.getVerbatim()
			print "Absolute difference (g):",abs(reactantsMass-productsMass)
			print "Mass reactants:"
			for RR in self.reactants:
				print RR.name, RR.mass
			print "Mass products:"
			for PP in self.products:
				print PP.name, PP.mass
			sys.exit()

	#********************
	def getVerbatim(self):
		if(self.verbatim!=None): return self.verbatim
		reactantsName = [x.name for x in self.reactants]
		productsName = [x.name for x in self.products]
		self.verbatim = (" + ".join(reactantsName))+" -> "+(" + ".join(productsName))
		return self.verbatim

	#********************
	def getVerbatimLatex(self):
		if(self.verbatimLatex!=None): return self.verbatimLatex
		reactantsName = [x.nameLatex for x in self.reactants]
		productsName = [x.nameLatex for x in self.products]
		self.verbatimLatex = (" + ".join(reactantsName))+" $\\to$ "+(" + ".join(productsName))
		return self.verbatimLatex

	#********************
	def getVerbatimHtml(self):
		if(self.verbatimHtml!=None): return self.verbatimHtml
		reactantsName = [x.nameHtml for x in self.reactants]
		productsName = [x.nameHtml for x in self.products]
		self.verbatimHtml = (" + ".join(reactantsName))+" &rarr; "+(" + ".join(productsName))
		return self.verbatimHtml

	#********************
	def getReactionHash(self):
		if(self.reactionHash!=None): return self.reactionHash
		reactantsName = sorted([x.nameFile for x in self.reactants])
		productsName = sorted([x.nameFile for x in self.products])
		self.reactionHash = ("_".join(reactantsName))+"__"+("_".join(productsName))
		return self.reactionHash

	#********************
	def getReactionHtmlRow(self,mySpecies=None):
		reactantsName = []
		for species in self.reactants:
			xspec = species.nameHref
			if(mySpecies!=None):
				if(species.name==mySpecies.name):
					xspec = "<b>"+species.nameHtml+"</b>"
			reactantsName.append(xspec)

		productsName = []
		for species in self.products:
			xspec = species.nameHref
			if(mySpecies!=None):
				if(species.name==mySpecies.name):
					xspec = "<b>"+species.nameHtml+"</b>"
			productsName.append(xspec)

		#reactantsName = sorted(reactantsName)
		#productsName = sorted(productsName)

		rpart = ("<td>+<td>".join(reactantsName))
		ppart = ("<td>+<td>".join(productsName))
		rspace = ("<td><td>".join([""]*(6-len(reactantsName))))
		pspace = ("<td><td>".join([""]*(10-len(productsName))))
		self.reactionHtmlRow = "<td>"+rpart+"<td>"+rspace+"<td>&rarr;<td>"+ppart+"<td>"+pspace
		self.reactionHtmlRow += "<td><a href=\"rate_"+self.getReactionHash()+".html\">details</a>"
		return self.reactionHtmlRow

	#********************
	def merge(self,myReaction):
		self.Tmin += myReaction.Tmin
		self.Tmax += myReaction.Tmax
		self.rate += myReaction.rate

	#*******************
	def plotRate(self,shortcuts,varRanges):
		self.evalRate(shortcuts,varRanges)
		self.doPlot()

	#********************
	#evaluate rates
	def evalRate(self,shortcuts,varRanges):

		#number of points
		imax = 100

		self.evaluation = []
		self.evalRate = []
		self.warnings = []
		self.shortcutsFound = dict()

		#loop on rates parts
		for icount in range(len(self.rate)):
			#get urrent rate
			rate = self.rate[icount]
			evaluation = dict()

			substFound = True
			#loop until shortcut found and replaced
			while(substFound):
				#remove spaces
				rate = rate.replace(" ","").lower()

				ops = ["+","-","/","*","(",")"]
				#add #s around operators
				for op in ops:
					rate = "#"+rate.replace(op,"#"+op+"#")+"#"

				#remove double exponent operator
				rate = rate.replace("#dexp#","#exp#")

				substFound = False
				#split rate at #s
				splitRate = [x for x in rate.split("#") if x!=""]
				#sort shortcuts by size
				var = sorted(shortcuts,key=lambda x:len(x),reverse=True)
				#loop on shortcuts to replace
				for variable in var:
					#loop on rate parts
					for i in range(len(splitRate)):
						#if shortcut found replace with expression
						if(splitRate[i]==variable.lower()):
							splitRate[i] = "("+shortcuts[variable]+")"
							self.shortcutsFound[variable] = shortcuts[variable]
							substFound = True
				#join rate back
				rate = ("".join(splitRate))

			#replace f90 numbers
			rate = rate.replace("d","e")

			ops = ["+","-","/","*","(",")"]
			for op in ops:
				rate = "#"+rate.replace(op,"#"+op+"#")+"#"

			#loop on available ranges
			for (variable,vrange) in varRanges.iteritems():
				#check if Tgas
				isTgas = (variable.lower()=="tgas")
				#get range limits
				(varMin,varMax) = vrange
				#log limits
				logVarMin = log10(varMin)
				logVarMax = log10(varMax)
				#create range
				vals = [i*(logVarMax-logVarMin)/(imax-1)+logVarMin for i in range(imax)]
				vals = [1e1**x for x in vals]

				#when Tgas check add limits if any
				if(isTgas):
					Tmin = varMin
					Tmax = varMax
					if(self.Tmin[icount]!=None): Tmin = float(utils.replaceTlims(self.Tmin[icount]))
					if(self.Tmax[icount]!=None): Tmax = float(utils.replaceTlims(self.Tmax[icount]))
					valsRange = [x for x in vals if(Tmin<=x and x<=Tmax)]
					valsRange = [Tmin]+valsRange+[Tmax]

				#when is not temperature if variable not in rate skip rate
				if(not(isTgas)):
					if(not("#"+variable.lower()+"#" in rate)): continue

				#store evaluated rate
				self.evalRate.append(rate)
				#store xdata and init ydata
				evaluation[variable] = dict()
				evaluation[variable]["xdata"] = vals
				evaluation[variable]["ydata"] = []
				#if Tgas store as interval
				if(isTgas):
					evaluation[variable]["xdataRange"] = valsRange
					evaluation[variable]["ydataRange"] = []
				hasEval = False
				#evaluate rate full range
				for val in vals:
					k = rate.replace("#"+variable.lower()+"#",str(val)).replace("#","")
					try:
						evaluation[variable]["ydata"].append(eval(k))
						hasEval = True
					except:
						evaluation[variable]["ydata"].append(None)

				#check if negative values found
				if(min(evaluation[variable]["ydata"])<0 and hasEval):
					if(isTgas):
						self.warnings.append("negative values when extrapolated")
					else:
						self.warnings.append("negative values found")

				#evaluate rate limited range
				if(isTgas and hasEval):
					kmin = eval(rate.replace("#"+variable.lower()+"#",str(Tmin)).replace("#",""))
					kmax = eval(rate.replace("#"+variable.lower()+"#",str(Tmax)).replace("#",""))

					evaluation[variable]["xlimits"] = [Tmin,Tmax]
					evaluation[variable]["ylimits"] = [kmin,kmax]
					for val in valsRange:
						k = rate.replace("#"+variable.lower()+"#",str(val)).replace("#","")
						try:
							evaluation[variable]["ydataRange"].append(eval(k))
							hasEval = True
						except:
							evaluation[variable]["ydataRange"].append(None)

					#check if negative values found (when limited range)
					if(min(evaluation[variable]["ydataRange"])<0 and hasEval):
						self.warnings.append("negative values found")


				#if not evaluated put none
				if(not(hasEval)):
					#add warning if rate not evaluated
					self.warnings.append("no rate evaluation")
					evaluation[variable] = None

			#store as class attribute
			self.evaluation.append(evaluation)

	#**************************
	def doPlot(self):
		import matplotlib.pyplot as plt
		#max orders of magnitude y axis
		yspanMax = 1e-6
		plt.clf()
		hasPlot = False
		#loop on different limited ranges
		for evaluation in self.evaluation:
			#get loop variable and evaluated rate
			for (variable,data) in evaluation.iteritems():
				if(data==None): continue
				hasPlot = True
				xdata = data["xdata"]
				ydata = data["ydata"]
				#plot full range
				plt.loglog(xdata,ydata,"r--")

				#if Tgas use limited ranges and plot limit points
				if(variable.lower()=="tgas"):
					xdataRange = data["xdataRange"]
					ydataRange = data["ydataRange"]
					plt.loglog(evaluation[variable]["xlimits"], evaluation[variable]["ylimits"],"ro")
					plt.loglog(xdataRange,ydataRange,"b-")
				else:
					xdataRange = xdata
					ydataRange = ydata

				#plot limited range
				plt.xlabel(variable)
				plt.ylabel("rate")
				plt.title(self.getVerbatimLatex())
				#set limits including max span
				plt.ylim(max(max(ydataRange)*yspanMax,min(ydataRange)*1e-1), max(ydataRange)*1e1)
				#set limits if constant
				if(min(ydataRange)==max(ydataRange)): plt.ylim(max(ydataRange)*1e-1,max(ydataRange)*1e1)

		#if value found save plot to png file
		if(hasPlot): plt.savefig("pngs/rate_"+str(self.getReactionHash())+"_"+variable+".png")


	#****************
	def makeHtmlPage(self,myOptions):

		fname = "htmls/rate_"+str(self.getReactionHash())+".html"

		reactantsList = [x.getHtmlName() for x in self.reactants]
		productsList = [x.getHtmlName() for x in self.products]

		table = []
		#table.append(["reactants", (", ".join(reactantsList))])
		#table.append(["products", (", ".join(productsList))])
		for icount in range(len(self.rate)):
			table.append(["rate", self.rate[icount]])
			table.append(["Tmin", self.Tmin[icount]])
			table.append(["Tmax", self.Tmax[icount]])

		for (shortcutName,shortcutExpression) in self.shortcutsFound.iteritems():
			table.append(["", shortcutName+" = "+shortcutExpression])

		for warning in self.warnings:
			table.append(["", warning])

		fout = open(fname,"w")
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">"+self.getVerbatimHtml()+"</p>\n")
		fout.write("<br>\n")
		fout.write("<a href=\"indexReactions.html\">back</a>\n")
		fout.write("<br><br>\n")
		fout.write("<table>\n")
		fout.write("<tr><th><th><th>\n")
		for (label,value) in table:
			if(value==None): continue
			separator = "&nbsp;&nbsp;:&nbsp;&nbsp;"
			if(label==""): separator = ""
			fout.write("<tr><td>"+label+"<td>"+separator+"<td>"+value+"\n")
		fout.write("<tr><th><th><th>\n")
		fout.write("</table>\n")

		for rng in myOptions.range:
			(rangeName,rangeValue) = [x.strip() for x in rng.split("=")]
			plotFileName = "pngs/rate_"+str(self.getReactionHash())+"_"+rangeName+".png"
			if(not(os.path.isfile(plotFileName))): continue
			fout.write("<img width=\"700px\" src=\"../"+plotFileName+"\">\n")

		fout.write(utils.getFile("footer.php"))
		fout.close()
