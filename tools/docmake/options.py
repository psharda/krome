class options:

	#constructor
	def __init__(self,fileName):

		optionsData = dict()

		optionsData["range"] = []

		#open options file
		fh = open(fileName,"rb")
		#loop on file lines
		for row in fh:
			srow = row.strip()
			#skip blank lines and comments
			if(srow==""): continue
			if(srow.startswith("#")): continue
			#get column-separated options
			(name,value) = [x.strip() for x in srow.split(":") if x!=""]

			if(name=="range"):
				optionsData[name].append(value)
			else:
				optionsData[name] = value

		for (k,v) in optionsData.iteritems():
			setattr(self,k,v)



