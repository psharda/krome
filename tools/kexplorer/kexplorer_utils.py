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

#********************
#Generate parenthesized contents in string as pairs (level, contents)
def parentheticContents(string):
	stack = []
	for i, c in enumerate(string):
	    if c == '{':
	        stack.append(i)
	    elif c == '}' and stack:
	        start = stack.pop()
	        yield (len(stack), string[start + 1: i])

#********************
#get content in parentses
def getParentheticContents(string):
	return list(parentheticContents(string))
