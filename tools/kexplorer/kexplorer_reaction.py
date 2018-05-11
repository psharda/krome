#kexplorer reaction class
class reaction:


	#*******************************
	#object constructor use reaction verbatim to init reactants and products
	def __init__(self,verbatim):


		#get reactants and products
		if verbatim.endswith("freeze-out"):
			RR = verbatim.split(" ")[0]
			PP = RR + "_grain"
		elif verbatim.endswith("evaporation"):
			PP = verbatim.split(" ")[0]
			RR = PP + "_grain"
		else:
			(RR,PP) = verbatim.split("->")

		self.reactants = [x for x in RR.split(" ") if(x!="" and x!="+")]
		self.products = [x for x in PP.split(" ") if(x!="" and x!="+")]

		self.reactantsHtml = [self.toHtml(x) for x in self.reactants]
		self.productsHtml = [self.toHtml(x) for x in self.products]

		#create DOT code for the given reaction
		self.dotCode = ""
		for reactant in self.reactants:
			for product in self.products:
				self.dotCode += "\""+reactant+"\" -> \""+product+"\";\n"

		#store verbatim
		self.verbatim = verbatim
		self.makeHtmlVerbatim()

		#init reaction data
		self.xvarData = []
		self.tgasData = []
		self.fluxData = []
		self.fluxNomrmMaxData = []
		self.fluxNormTotData = []

	#*******************
	#append flux data
	def addData(self, xvar, Tgas, flux, fluxNormMax, fluxNormTot):
		self.xvarData.append(xvar)
		self.tgasData.append(Tgas)
		self.fluxData.append(flux)
		self.fluxNomrmMaxData.append(fluxNormMax)
		self.fluxNormTotData.append(fluxNormTot)


	#*****************
	def toHtml(self,species):
		repSup = ["+","-"]
		for rep in repSup:
			species = species.replace(rep,"<sup>"+rep+"</sup>")
		aspec = list(species)

		species = ""
		for x in aspec:
			try:
				float(x)
				species += "<sub>"+x+"</sub>"
			except:
				species += x

		return species

	#*****************
	def makeHtmlVerbatim(self):

		htmlRR = ("<td>+&nbsp;<td>".join(self.reactantsHtml))
		htmlPP = ("<td>+&nbsp;<td>".join(self.productsHtml))

		tdmax = 4
		htmlRR += "<td>"*(tdmax-htmlRR.count("<td>"))
		htmlPP += "<td>"*(tdmax-htmlPP.count("<td>"))

		htmlRR = "<td>"+htmlRR
		htmlPP = "<td>"+htmlPP

		self.htmlVerbatim = htmlRR+"<td>&rarr;"+htmlPP


