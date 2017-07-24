import glob

#folder name
folder = "networks/"

#parse format reactants and products
#returns indexes with positions
def parseFormat(arow):
	idxs = []
	for i in range(len(arow)):
		if(arow[i] in ["R","P"]): idxs.append(i)
	return idxs

#check if variable is an integer number
def isNumber(arg):
	try:
		int(arg)
		return True
	except:
		return False

#default format list
defaultFormat = ("idx,R,R,R,P,P,P,P").split(",")

#loop on files in folder
for ntw in sorted(glob.glob(folder+"*")):
	idxs = parseFormat(defaultFormat)

	print "****************"
	print ntw
	reactionCount = 0
	species = []
	#loop on file
	for row in open(ntw,"rb"):
		srow = row.strip()
		#skip comments and blanks
		if(srow==""): continue
		if(srow.startswith("#")): continue
		if(srow.startswith("!")): continue
		arow = srow.split(",")
		#get format
		if(srow.startswith("@format:")):
			idxs = parseFormat(arow)
		#skip lines not starting with integers (index)
		if(not(isNumber(arow[0]))): continue

		#loop on indexes with reactants/products
		for idx in idxs:
			#print idx,arow
			species.append(arow[idx])
		reactionCount += 1

	#sort and unique species
	species = sorted(list(set(species)))
	species = [x for x in species if(x!="")]

	#print output
	print "Number of reactions:",reactionCount
	print "Number of species:",len(species)
	print (", ".join(species))
	print


