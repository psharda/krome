#kexplorer reaction class
class reaction:


	#*******************************
	#object constructor use reaction verbatim to init reactants and products
	def __init__(self,verbatim):

		#get reactants and products
		(RR,PP) = verbatim.split("->")
		self.reactants = [x for x in RR.split(" ") if(x!="" and x!="+")]
		self.products = [x for x in PP.split(" ") if(x!="" and x!="+")]

		#create DOT code for the given reaction
		self.dotCode = ""
		for reactant in self.reactants:
			for product in self.products:
				self.dotCode += "\""+reactant+"\" -> \""+product+"\";\n"

		#store verbatim
		self.verbatim = verbatim

		#init reaction data
		self.xvarData = []
		self.fluxData = []
		self.fluxNomrmMaxData = []
		self.fluxNormTotData = []

	#*******************
	#append flux data
	def addData(self, xvar, flux, fluxNormMax, fluxNormTot):
		self.xvarData.append(xvar)
		self.fluxData.append(flux)
		self.fluxNomrmMaxData.append(fluxNormMax)
		self.fluxNormTotData.append(fluxNormTot)

