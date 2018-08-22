#kexplorer element class
class element:


	#*******************************
	#object constructor
	def __init__(self):

		#init reaction data
		self.xvarData = []
		self.tgasData = []
		self.timeData = []
		self.abundanceData = []
		self.blockCnt = -1


	#*******************
	#append time evolution of element
	def addData(self, time, Tgas, xvar, abundance, newBlock):
		#append list to list for new block in input file
		if(newBlock):
			self.timeData.append([])
			self.tgasData.append([])
			self.xvarData.append([])
			self.abundanceData.append([])
			self.blockCnt += 1

		idx = self.blockCnt
		self.timeData[idx].append(time)
		self.tgasData[idx].append(Tgas)
		self.xvarData[idx].append(xvar)
		self.abundanceData[idx].append(abundance)

	#*******************
	# sort on tgas and xvar
	def sortData(self):
		zipped = zip(self.tgasData, self.xvarData, self.timeData, self.abundanceData)
		sorted_zipped = sorted(zipped, key=lambda x: (x[1], x[0]))
		(self.tgasData,
		self.xvarData,
		self.timeData,
		self.abundanceData) = (list(t) for t in zip(*sorted_zipped))
