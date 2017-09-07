#useful functions

#*********************
#check if arg is a number
def isNumber(arg):
	try:
		float(arg)
		return True
	except ValueError:
		return False

#********************
#get shortcuts of temperature that can be used in network
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

#********************
#uses empty line as separator
def blockSeparator(line):
	return line=='\n'
