import sys,species,utils,os,urllib
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
	#parse csv reaction file row (constructor)
	def __init__(self,row,reactionFormat,atomSet,reactionType):
		if(not(reactionFormat.startswith("@format:"))):
			sys.exit("ERROR: wrong format "+reactionFormat)
		splitFormat = reactionFormat.replace("@format:","").split(",")
		splitFormat = [x.lower() for x in splitFormat]
		arow = row.strip().split(",",len(splitFormat))
		arow = [x.strip() for x in arow]

		specials = ["","G"]

		self.rate = []
		self.reactants = []
		self.products = []
		self.Tmin = [None]
		self.Tmax = [None]
		self.reactionType = reactionType

		#loop on format parts (and parse species)
		for i in range(len(splitFormat)):
			part = splitFormat[i]
			#reaction index
			if(part=="idx"):
				self.index = int(arow[i])
			#reactant
			elif(part=="r"):
				if(arow[i].upper() in specials): continue
				spec = species.species(arow[i],atomSet)
				self.reactants.append(spec)
			#product
			elif(part=="p"):
				if(arow[i].upper() in specials): continue
				spec = species.species(arow[i],atomSet)
				self.products.append(spec)
			#min temperature
			elif(part=="tmin"):
				self.Tmin = [arow[i].replace("d","e")]
				if(arow[i].lower()=="none" or arow[i]==""): self.Tmin = [None]
			#max temperature
			elif(part=="tmax"):
				self.Tmax = [arow[i].replace("d","e")]
				if(arow[i].lower()=="none" or arow[i]==""): self.Tmax = [None]
			#rate coefficient
			elif(part=="rate"):
				self.rate.append(arow[i])
			else:
				print "ERROR: unknow format element "+part
				sys.exit(reactionFormat)

		#add cosmic ray if not present
		hasCR = ("CR" in [x.name for x in self.reactants])
		if(not(hasCR) and reactionType=="CR"):
				spec = species.species("CR",atomSet)
				self.reactants.append(spec)

		#check mass and charge conservation
		self.check()

	#**************
	#check reaction charge and mass conservation
	def check(self):
		#check charge
		reactantsCharge = sum([x.charge for x in self.reactants])
		productsCharge = sum([x.charge for x in self.products])
		if(reactantsCharge!=productsCharge):
			print "ERROR: charge problems"
			print self.getVerbatim()
			sys.exit()

		#check mass
		reactantsMass = sum([x.mass for x in self.reactants])
		productsMass = sum([x.mass for x in self.products])
		me = 9.10938356e-28 #electron mass, g
		#error if difference in mass is greater than half electron mass
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
	#get verbatim reaction, e.g.H2+H2->H+H+H2
	def getVerbatim(self):
		if(self.verbatim!=None): return self.verbatim
		reactantsName = [x.name for x in self.reactants]
		productsName = [x.name for x in self.products]
		self.verbatim = (" + ".join(reactantsName))+" -> "+(" + ".join(productsName))
		return self.verbatim

	#********************
	#get reaction verbatim in latex format, e.g. H$^+$+e$^-$$\to$H
	def getVerbatimLatex(self):
		if(self.verbatimLatex!=None): return self.verbatimLatex
		reactantsName = [x.nameLatex for x in self.reactants]
		productsName = [x.nameLatex for x in self.products]
		self.verbatimLatex = (" + ".join(reactantsName))+" $\\to$ "+(" + ".join(productsName))
		return self.verbatimLatex

	#********************
	#get HTML verbatim, e.g. H<sub>2</sub>&rarr;H+H
	def getVerbatimHtml(self):
		if(self.verbatimHtml!=None): return self.verbatimHtml
		reactantsName = [x.nameHtml for x in self.reactants]
		productsName = [x.nameHtml for x in self.products]
		self.verbatimHtml = (" + ".join(reactantsName))+" &rarr; "+(" + ".join(productsName))
		return self.verbatimHtml

	#********************
	#get reaction unique hash, e.g. H_H__H2
	def getReactionHash(self):
		if(self.reactionHash!=None): return self.reactionHash
		reactantsName = sorted([x.nameFile for x in self.reactants])
		productsName = sorted([x.nameFile for x in self.products])
		self.reactionHash = ("_".join(reactantsName))+"__"+("_".join(productsName))
		return self.reactionHash

	#********************
	#get html table row with bold mySpecies when present
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

		rpart = ("<td>+<td>".join(reactantsName))
		ppart = ("<td>+<td>".join(productsName))
		#fill empty spaces
		rspace = ("<td><td>".join([""]*(6-len(reactantsName))))
		pspace = ("<td><td>".join([""]*(10-len(productsName))))

		self.reactionHtmlRow = "<td>&nbsp;"+rpart+"<td>"+rspace+"<td>&rarr;<td>"+ppart+"<td>"+pspace
		self.reactionHtmlRow += "<td><a href=\"rate_"+self.getReactionHash()+".html\">details</a>"

		return self.reactionHtmlRow


	#*******************
	def getAtoms(self):
		atoms = [x for x in self.getSpecies()]
		return list(set(atoms))

	#********************
	#get list of species
	def getSpecies(self):
		return self.reactants+self.products

	#********************
	#merge limits and rates with another rate
	def merge(self,myReaction):
		self.Tmin += myReaction.Tmin
		self.Tmax += myReaction.Tmax
		self.rate += myReaction.rate

	#*******************
	#plot rate coefficient
	def plotRate(self,shortcuts,varRanges):
		self.evalRate(shortcuts,varRanges)
		self.doPlot()
		self.saveEvals()

	#********************
	#search for rate variables in the rate
	def getRateVariables(self,myOptions):
		rateVariables = []
		#loop on user defined variables
		for variable in myOptions.getRanges().keys():
			#loop on rates to search variable name
			for rate in self.rate:
				#append when variable found
				if(variable in rate.lower()): rateVariables.append(variable)
		#if no variable found assumes Tgas
		if(len(rateVariables)==0): rateVariables = ["tgas"]

		#return unique list
		return list(set(rateVariables))

	#*****************
	#boolean returned if variable found in rates
	def hasVariable(self,myOptions,variable):
		return (variable.lower() in [x.lower() for x in self.getRateVariables(myOptions)])

	#********************
	#evaluate rates
	def evalRate(self,shortcuts,varRanges):

		#number of points in the plot
		imax = 100

		self.evaluation = []
		self.evalRate = []
		self.warnings = []
		self.shortcutsFound = dict()

		#F90 expected operators in F90 rate expression
		ops = ["+","-","/","*","(",")"]


		#loop on rates parts
		for icount in range(len(self.rate)):
			#get current rate
			rate = self.rate[icount]
			evaluation = dict()

			substFound = True
			#loop until shortcut found and replaced
			while(substFound):
				#remove spaces
				rate = rate.replace(" ","").lower()

				#surround F90 operators with # symbols
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

			#replace F90 numbers for evaluation, d->e
			rate = rate.replace("d","e")

			#surround F90 operators with # symbols
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
				#create variable range
				vals = [i*(logVarMax-logVarMin)/(imax-1)+logVarMin for i in range(imax)]
				vals = [1e1**x for x in vals]

				#when Tgas check add limits if any
				if(isTgas):
					Tmin = varMin
					Tmax = varMax
					if(self.Tmin[icount]!=None): Tmin = float(utils.replaceTlims(self.Tmin[icount]))
					if(self.Tmax[icount]!=None): Tmax = float(utils.replaceTlims(self.Tmax[icount]))
					Tmin = max(varMin,Tmin)
					Tmax = min(varMax,Tmax)
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

				#evaluation exist flag
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
	#do plot (PNG)
	def doPlot(self):
		import matplotlib.pyplot as plt

		plt.clf()
		#max orders of magnitude y axis
		yspanMax = 1e-10
		hasPlot = False
		#loop on different limited ranges
		for evaluation in self.evaluation:
			#get loop variable and evaluated rate
			for (variable,data) in evaluation.iteritems():

				if(data==None): continue
				hasPlot = True
				xdata = data["xdata"]
				ydata = data["ydata"]

				if(variable.lower()=="tgas"):
					#plot full range
					plt.loglog(xdata,ydata,"r--")

					#if Tgas use limited ranges and plot limit points
					xdataRange = data["xdataRange"]
					ydataRange = data["ydataRange"]
					plt.loglog(evaluation[variable]["xlimits"], evaluation[variable]["ylimits"],"ro")
					plt.loglog(xdataRange,ydataRange,"b-")
				else:
					plt.clf()
					xdataRange = xdata
					ydataRange = ydata
					plt.loglog(xdataRange,ydataRange)


				plt.grid(b=True, color='0.65',linestyle='--')
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
	#save evaluation as a json structure
	def saveEvals(self):
		import json

		#json file name
		fname = "evals/rate_"+str(self.getReactionHash())+".json"

		#convert to json
		jsonData = json.dumps(self.evaluation)

		#dump to file
		fout = open(fname,"w")
		fout.write(jsonData)
		fout.close()


	#****************
	#make corresponding HTML page
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
		fout.write("<br>\n")
		urlencoded = urllib.quote_plus(" + ".join([x.name for x in self.reactants]))
		urlkida = "http://kida.obs.u-bordeaux1.fr/search.html?species="+urlencoded+"&reactprod=both&astroplaneto=Both&ionneutral=ion&isomers=1&ids="
		fout.write("<a href=\""+urlkida+"\" target=\"_blank\">search in KIDA</a>\n")
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

		fout.write(utils.getFooter("footer.php"))
		fout.close()

