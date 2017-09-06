#useful functions

#*********************
def isNumber(arg):
	try:
		float(arg)
		return True
	except ValueError:
		return False

#********************
def getShortcuts():
	shortcut = []
	fileName = "temperatureShortcuts.dat"
	fh = open(fileName,"rb")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		(variable,expression) = [x.strip() for x in srow.split("!")[0].split("=")]
		shortcut.append((variable,expression))
	return shortcut
