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

#********************
#find all indices of an element in a list
def indicesElemList(list,elem):
	indices = [i for i, x in enumerate(list) if x == elem]
	if indices == []:
		indices = [None]

	return indices

#********************
#sub and super script to LaTeX format
def subSuper2latex(name):
	#maks subcripts of numbers
	for part in name:
		if(isNumber(part)):
			name = name.replace(part,"$_{"+part+"}$")

	#make superscript of all ions
	for sign in ["+","-"]:
		if sign in name:
			allSigns = name.count(sign)*sign
			name = name.replace(allSigns,"$^{"+allSigns+"}$")
	return name

#********************
#specials to LaTeX format
def special2latex(name):
	if name =="E": name = "e$^-$"
	if name =="g": name = "$\\gamma$"
	if "_c" in name: name = name.replace("_c","$_c$")
	if "_dust" in name: name = name.replace("_dust","$_{dust}$")
	return name

#********************
#limits to LaTeX format
def limits2latex(name):
	name = name.replace("<=", " $\leqslant$ ").replace(">="," $\geqslant$ ")
	name = name.replace("<","$ < $").replace(">"," $>$ ")
	return name

#********************
#simplify KROME limits
def limitsSimple(name):
	name = name.replace(".LE.","<=").replace(".GE.",">=")
	name = name.replace(".LT.","<").replace(".GT.",">")
	name = name.replace("NONE","")
	name = name.replace("d","e")
	return name

#********************
#character to float
def char2float(arg):
	if isNumber(arg):
		return float(arg)
	else:
		return arg

#********************
#character to int
def char2int(arg):
	if isNumber(arg):
		return int(float(arg))
	else:
		return arg
