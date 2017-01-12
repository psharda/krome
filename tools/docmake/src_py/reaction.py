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
	safeExtrapolate = dict()

	#********************
	#parse csv reaction file row (constructor)
	def __init__(self,row,reactionFormat,atomSet,reactionType):

		if(reactionType=="KIDA"):
			self.parseFormatKIDA(row,reactionFormat,atomSet,reactionType)
		else:
			self.parseFormatKROME(row,reactionFormat,atomSet,reactionType)

	#********************
	def parseFormatKROME(self,row,reactionFormat,atomSet,reactionType):

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

		#add cosmic rays if not present
		hasCR = ("CR" in [x.name for x in self.reactants])
		if(not(hasCR) and reactionType=="CR"):
				spec = species.species("CR",atomSet)
				self.reactants.append(spec)

		#check mass and charge conservation
		self.check()


	#********************
	def parseFormatKIDA(self,row,reactionFormat,atomSet,reactionType):


		specials = ["","G","PHOTON"]

		#variables for cosmic rays and photochemistry
		CRvar = "user_crate" #name of the CR flux variable
		Avvar = "user_Av" #name of the Av variable

		#number of reactants and products expected in KIDA file
		maxReactants = 3
		maxProducts = 5

		#spacing format
		fmt = [11]*3 + [1] + 5*[11] +[1] + 3*[11] + [8,9] + [1,4,3] + 2*[7] + [3,6,3,2]
		#keys names
		keys = ["R"+str(i) for i in range(maxReactants)] + ["x"] \
			+ ["P"+str(i) for i in range(maxProducts)] +["x"] \
			+ ["a","b","c"] + ["F","g"] + ["x","unc","type"] \
			+ ["tmin","tmax"] + ["formula","num","subnum","recom"]

		self.rate = []
		self.reactants = []
		self.products = []
		self.Tmin = [None]
		self.Tmax = [None]
		self.reactionType = reactionType

		srow = row.strip()

		dataRow = dict()
		position = 0
		#loop on format to get data from the row as a dictionary
		for i in range(len(fmt)):
			dataRow[keys[i]] = srow[position:position+fmt[i]].strip()
			position += fmt[i]

		self.Tmin = [dataRow["tmin"]]
		self.Tmax = [dataRow["tmax"]]

		for i in range(maxReactants):
			v = dataRow["R"+str(i)].strip()
			if(v=="e-"): v = "E"
			if(v.upper() in specials): continue
			spec = species.species(v,atomSet)
			self.reactants.append(spec)

		for i in range(maxProducts):
			v = dataRow["P"+str(i)].strip()
			if(v=="e-"): v = "E"
			if(v.upper() in specials): continue
			spec = species.species(v,atomSet)
			self.products.append(spec)

		#Formula is a number that referes to the formula needed to compute the rate coefficient of the reaction.
		#see http://kida.obs.u-bordeaux1.fr/help
		#1: Cosmic-ray ionization (direct and undirect)
		#2: Photo-dissociation (Draine)
		#3: Kooij
		#4: ionpol1
		#5: ionpol2
		#6: Troe fall-off (NOT SUPPORTED!)
		arow = dataRow
		arow["formula"] = int(arow["formula"])
		if(arow["formula"]==0):
			KK = "auto"
		elif(arow["formula"]==1):
			KK = arow["a"]+"*"+CRvar
		elif(arow["formula"]==2):
			KK = arow["a"]
			if(float(arow["c"])!=0e0): KK += "*exp(-"+arow["c"]+"*"+Avvar+")"
		elif(arow["formula"]==3):
			KK = arow["a"]
			if(float(arow["b"])!=0e0): KK += "*(T32)**("+arow["b"]+")"
			if(float(arow["c"])!=0e0): KK += "*exp(-"+arow["c"]+"*invT)"
		elif(arow["formula"]==4):
			KK = arow["a"]
			if(float(arow["b"])!=1e0): KK += "*"+arow["b"]
			gpart = ""
			if(float(arow["c"])!=0e0): gpart = "+ 0.4767d0*"+arow["c"]+"*sqrt(3d2*invT)"
			KK += "*(0.62d0 "+gpart+")"
		elif(arow["formula"]==5):
			KK = arow["a"]
			if(float(arow["b"])!=1e0): KK += "*"+arow["b"]
			gpart = ""
			if(float(arow["c"])!=0e0):
				gpart = "+ 0.0967d0*"+arow["c"]+"*sqrt(3d2*invT) + "
				gpart += arow["c"]+"**2*28.501d0*invT"
				KK += "*(1d0 "+gpart+")"
		else:
			print "ERROR: KIDA formula "+str(arow["formula"])+" not supported!"

		KK = KK.replace("--","+").replace("++","+").replace("-+","-").replace("+-","-")

		self.rate.append(KK)


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
	#get reaction NON-unique (unsorted) hash, e.g. H_H__H2
	def getReactionHashUnsorted(self):
		reactantsName = [x.nameFile for x in self.reactants]
		productsName = [x.nameFile for x in self.products]
		return ("_".join(reactantsName))+"__"+("_".join(productsName))

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
		self.evaluateExtrapolation(varRanges)
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
		import matplotlib
		#try to load AGG for PNG rendering (slightly faster)
		try:
			matplotlib.use('AGG')
		except:
			pass

		import matplotlib.pyplot as plt

		#turn off interactivity
		plt.ioff()

		#cancel current plot
		plt.clf()
		#max orders of magnitude y axis
		yspanMax = 1e-10
		hasPlot = False
		ydataAll = []
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
					ydataAll += ydataRange
					plt.loglog(evaluation[variable]["xlimits"], evaluation[variable]["ylimits"],"ro")
					plt.loglog(xdataRange,ydataRange,"b-")
				else:
					plt.clf()
					xdataRange = xdata
					ydataRange = ydata
					ydataAll += ydataRange
					plt.loglog(xdataRange,ydataRange)

		#plot only if data are available
		if(hasPlot):
			plt.grid(b=True, color='0.65',linestyle='--')
			#plot limited range
			plt.xlabel(variable)
			plt.ylabel("rate")
			plt.title(self.getVerbatimLatex())
			#set limits including max span
			plt.ylim(max(max(ydataAll)*yspanMax,min(ydataAll)*1e-1), max(ydataAll)*1e1)
			#set limits if constant
			if(min(ydataAll)==max(ydataAll)): plt.ylim(max(ydataAll)*1e-1,max(ydataAll)*1e1)

			#if value found save plot to png file
			plt.savefig("pngs/rate_"+str(self.getReactionHash())+"_"+variable+".png", dpi=150)

	#******************
	#evaluate rate extrapolation for the current reaction
	def evaluateExtrapolation(self,varRanges):

		#init Tgas limits
		xMin = 1e99
		xMax = -1e99
		#init flags
		hasData = isIncreasingMin = isDecreasingMax = False
		isAlwaysPositiveMin = isAlwaysPositiveMax = False
		#loop on different limited ranges
		for evaluation in self.evaluation:
			#get only Tgas data
			for (variable,vdata) in evaluation.iteritems():
				if(variable.lower()!="tgas"): continue
				#store data
				data = vdata
				#store ranges
				varRange = varRanges[variable]

			#skip missing data
			if(data==None): continue
			hasData = True
			#copy data locally
			xdataRange = data["xdataRange"]
			ydataRange = data["ydataRange"]
			xdata = data["xdata"]
			ydata = data["ydata"]
			#number of data points
			ndata = len(xdata)

			#check smaller rate interval
			if(min(xdataRange)<xMin):
				#store min value
				xMin = min(xdataRange)
				#get ydata outside interval
				ydataOutside = [ydata[i] for i in range(ndata) if(xdata[i]<xMin)]
				#check if ydata are increasing outside
				isIncreasingMin = (sorted(ydataOutside)==ydataOutside)
				#check if data are always positive outside (only if data are present)
				isAlwaysPositiveMin = True
				if(len(ydataOutside)>0): isAlwaysPositiveMin = (min(ydataOutside)>0e0)
				#store min Tgas in data structure
				self.safeExtrapolate["Tmin"] = xMin

			if(max(xdataRange)>xMax):
				xMax = max(xdataRange)
				ydataOutside = [ydata[i] for i in range(ndata) if(xdata[i]>xMax)]
				isDecreasingMax = (sorted(ydataOutside)==ydataOutside[::-1])
				isAlwaysPositiveMax = True
				if(len(ydataOutside)>0): isAlwaysPositiveMax = (min(ydataOutside)>0e0)
				self.safeExtrapolate["Tmax"] = xMax

		#check extrapolation only if has data
		if(hasData):
			#store extrapolated Tgas limits
			self.safeExtrapolate["TminExtrapolated"] = min(varRange)
			self.safeExtrapolate["TmaxExtrapolated"] = max(varRange)
			#store if extrapolation is safe or not
			self.safeExtrapolate["lower"] = (isAlwaysPositiveMin and isIncreasingMin)
			self.safeExtrapolate["upper"] = (isAlwaysPositiveMax and isDecreasingMax)

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

		header = "<tr><th><th><th>\n"

		table = []
		tableNotes = []
		#table.append(["reactants", (", ".join(reactantsList))])
		#table.append(["products", (", ".join(productsList))])
		for icount in range(len(self.rate)):
			table.append(["header", header])
			table.append(["rate", self.rate[icount]])
			table.append(["Tmin", self.Tmin[icount]])
			table.append(["Tmax", self.Tmax[icount]])

		for (shortcutName,shortcutExpression) in self.shortcutsFound.iteritems():
			tableNotes.append(shortcutName+" = "+shortcutExpression)

		for warning in self.warnings:
			tableNotes.append(warning)

		fout = open(fname,"w")
		fout.write(utils.getFile("header.php"))
		fout.write("<p style=\"font-size:30px\">"+self.getVerbatimHtml()+"</p>\n")
		fout.write("<br>\n")
		fout.write("<a href=\"indexReactions.html\">back</a><br>\n")
		urlencoded = urllib.quote_plus(" + ".join([x.name for x in self.reactants]))
		urlkida = "http://kida.obs.u-bordeaux1.fr/search.html?species="+urlencoded+"&reactprod=reactants&astroplaneto=Both&ionneutral=ion&isomers=1&ids="
		fout.write("<a href=\""+urlkida+"\" target=\"_blank\">search in KIDA</a><br>\n")
		urlJSON = "../evals/rate_"+str(self.getReactionHash())+".json"
		fout.write("<a href=\""+urlJSON+"\">get rate evaluation in JSON format</a>\n")
		fout.write("<br><br>\n")
		fout.write("<table>\n")
		for (label,value) in table:
			if(value==None): continue
			if(label=="header"):
				fout.write(value)
			else:
				separator = "&nbsp;&nbsp;:&nbsp;&nbsp;"
				fout.write("<tr><td>"+label+"<td>"+separator+"<td>"+value+"\n")
		fout.write(header)
		fout.write("</table><br><br>\n")

		bulletPoint = "&nbsp;&nbsp;&#9656;&nbsp;"
		#add notes and warnings
		if(len(tableNotes)>0):
			fout.write("Notes and warnings:<br>\n")
			for value in tableNotes:
				if(value!=None): fout.write(bulletPoint+value+"<br>\n")


		#extrapolation lower limit
		TminExtrapolated = utils.htmlExpBig(self.safeExtrapolate["TminExtrapolated"])
		Tmin = utils.htmlExpBig(self.safeExtrapolate["Tmin"])
		extrapolCheckMin = ("SAFE" if self.safeExtrapolate["lower"] else "NOT SAFE")
		if(self.safeExtrapolate["TminExtrapolated"]!=self.safeExtrapolate["Tmin"]):
			fout.write(bulletPoint+"Extrapolation in range ["+TminExtrapolated+", "+Tmin+"] K is <b>"+extrapolCheckMin+"</b><br>")

		#extrapolation upper limit
		TmaxExtrapolated = utils.htmlExpBig(self.safeExtrapolate["TmaxExtrapolated"])
		Tmax = utils.htmlExpBig(self.safeExtrapolate["Tmax"])
		extrapolCheckMax = ("SAFE" if self.safeExtrapolate["upper"] else "NOT SAFE")
		if(self.safeExtrapolate["TmaxExtrapolated"]!=self.safeExtrapolate["Tmax"]):
			fout.write(bulletPoint+"Extrapolation in range ["+Tmax+", "+TmaxExtrapolated+"] K is <b>"+extrapolCheckMax+"</b><br>")

		for rng in myOptions.range:
			(rangeName,rangeValue) = [x.strip() for x in rng.split("=")]
			plotFileName = "pngs/rate_"+str(self.getReactionHash())+"_"+rangeName+".png"
			if(not(os.path.isfile(plotFileName))): continue
			fout.write("<img width=\"700px\" src=\"../"+plotFileName+"\">\n")

		fout.write(utils.getFooter("footer.php"))
		fout.close()

