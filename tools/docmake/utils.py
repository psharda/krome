#***********************
def getFile(fileName):
	with open(fileName, 'r') as content_file:
		return content_file.read()


#*********************
def getAtomSet(fileName):
	atomSet = dict()

	refMass = dict()
	refMass["me"] = 9.10938356e-28 #electron mass, g
	refMass["mp"] = 1.6726219e-24 #proton mass, g
	refMass["mn"] = 1.6749286e-24 #neutron mass, g
	refMass["mep"] = refMass["me"]+refMass["mp"]
	refMass["mepn"] = refMass["mep"]+refMass["mn"]

	fh = open(fileName,"rb")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		arow = srow.replace("\t"," ").split(" ")
		arow = [x.strip() for x in arow if x.strip()!=""]
		reps = sorted(refMass.keys(), key=lambda x:len(x),reverse=True)
		for rep in reps:
			arow[1] = arow[1].replace(rep,str(refMass[rep]))
		atomSet[arow[0].upper()] = eval(arow[1])
	return atomSet

#*********************
def replaceTlims(arg):
	reps = ["<",">",".LE.",".GE.",".LT.",".GT.","<=",">="]
	reps = sorted(reps,key=lambda x:len(x),reverse=True)
	for rep in reps:
		arg = arg.lower().replace(rep.lower(),"")
	return arg.replace("d","e")

#*********************
def isNumber(arg):
	try:
		float(arg)
		return True
	except ValueError:
		return False

#********************
def getShortcuts():
	shortcut = dict()
	fileName = "shortcuts.dat"
	fh = open(fileName,"rb")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		(variable,expression) = [x.strip() for x in srow.split("=")]
		shortcut[variable] = expression
	return shortcut

