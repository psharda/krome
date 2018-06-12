import os,datetime,math,re,regex

#***********************
def getHtmlProperty(item):
	if(item=="tableRowBgcolor"):
		return  "#D3D3D3"
	else:
		print "ERROR: item "+item+" unknown!"
		sys.exit()

#***********************
def getFile(fileName):
	with open(fileName, 'r') as content_file:
		return content_file.read()

#***********************
def getFooter(sourceFile):
	return getFile("footer.php").replace("#FOOTER_INFO#", getFooterInfo())

#**********************
def getChangeset():
	#name of the git master file
	masterfile = "../../.git/refs/heads/master"
	changeset = ("x"*7) #default unknown changeset
	#if git master file exists grep the changeset
	if(os.path.isfile(masterfile)):
		changeset = open(masterfile,"rb").read()
	return changeset.strip()

#***********************
def getCurrentTime():
	return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

#***********************
def getFooterInfo():

	changeset = getChangeset()
	bitbucketLink = "https://bitbucket.org/tgrassi/krome/commits/"+changeset
	bitbucketLanding = "<a href=\"https://bitbucket.org/tgrassi/krome\" target=\"_blank\">KROME</a>"
	hrefWiki = "<a href=\"https://bitbucket.org/tgrassi/krome/wiki/docmaker\" target=\"_blank\">docmake</a>"
	return "documentation generated with "+hrefWiki+" ("+bitbucketLanding+") - changeset: <a href=\""\
		+ bitbucketLink + "\" target=\"_blank\">" + changeset[:7] + "</a> - " + getCurrentTime()


#*********************
#load polarizability data from file
def getPolarizabilityData(fileName):

	print "loading polarizability data from "+fileName

	#data dictionary
	polData = dict()

	#open file to read
	fh = open(fileName,"rb")
	#loop on file
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		srow = srow.replace("\t"," ")
		(species,value) = [x.strip() for x in srow.split(" ") if(x!="")]
		polData[species] = float(value)*1e-24 #AA3->cm3

	fh.close()

	return polData


#*********************
#load thermochemical data from a burcat-like file
def getThermochemicalData(fileName):

	cnum = 5 #number of coefficients per line
	clen = 15 #number of characters per coefficient

	print "loading thermochemical data from "+fileName

	#thermochemical data dictionary
	thermochemDict = dict()

	#open file to read
	fh = open(fileName,"rb")
	#loop on file
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		arow = [x for x in srow.split(" ") if(x!="")]
		#last column is an integer (line count)
		lineIdx = int(arow[-1])
		#first line starts with species name
		if(lineIdx==1):
			#species name
			specName = arow[0]
			#temperature limits from first line
			(Tmin,Tmax,Tmid) = [float(x) for x in arow[-4:-1]]
			#init and store temperature limits in the dictionary
			thermochemDict[specName] = dict()
			thermochemDict[specName]["Tmin"] = Tmin
			thermochemDict[specName]["Tmid"] = Tmid
			thermochemDict[specName]["Tmax"] = Tmax
			coefs = []
		elif(lineIdx>1 and lineIdx<5):
			#loop on coefficients and store
			for i in range(cnum):
				coefs.append(float(row[i*clen:(i+1)*clen]))
			#if last line store coefficients
			if(lineIdx==4):
				thermochemDict[specName]["lowT"] = coefs[7:]
				thermochemDict[specName]["highT"] = coefs[:7]
		else:
			print "ERROR: unkonw reading index in "+fileName+" for line:"
			print srow
			sys.exit()
	fh.close()

	return thermochemDict


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
		atomSet[arow[0]] = eval(arow[1])

	fh.close()

	return atomSet

#*********************
#html format only if geq than maxlim
def htmlExpBig(arg,digits=2,maxLim=1e3):
	if(arg==None): return str(arg)
	if(arg==0): return "0"
	if(arg>=maxLim):
		return htmlExp(arg,digits=digits)
	else:
		return str(arg)

#*********************
def htmlExp(arg,digits=2):
	if(arg==None): return str(arg)
	if(arg==0): return "0"
	xp = int(math.log10(abs(arg)))-1
	mt = arg/1e1**(xp)
	mt = int(mt*1e1**digits)/1e1**digits
	if(mt==1e0):
		return "10<sup>"+str(xp)+"</sup>"
	if((mt!=1e0) and (mt==int(mt))):
		return str(int(mt))+"&times;10<sup>"+str(xp)+"</sup>"
	else:
		return str(mt)+"&times;10<sup>"+str(xp)+"</sup>"

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
#character to int
def char2int(arg, when_below=1e4):
	if isNumber(arg):
                f = float(arg)
                if(f < when_below):
		        return int(float(arg))
	return arg

#********************
def getShortcuts():
	shortcut = dict()
	fileName = "shortcuts.dat"
        absPath = os.path.join(os.path.dirname(__file__), "..", fileName)
        absPath = os.path.abspath(absPath)
	fh = open(absPath,"rb")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		(variable,expression) = [x.strip() for x in srow.split("=")]
		shortcut[variable] = expression
	return shortcut

#********************
#get shortcuts of temperature that can be used in network
def getShortcutsLatex():
        from options import latexoptions
        return latexoptions.substitutions

#********************
#check if variable is already a temperature shortcut
def isTemperatureShortcut(var):
	shortcuts = getShortcutsLatex()
	if var in [item[0] for item in shortcuts]:
		return True
	else:
		return False

#********************
#Generate parenthesized contents in string as pairs (level, contents)
def parentheticContents(string, parenteses):
	opened = parenteses[0]
	closed = parenteses[1]
	stack = []
	for i, c in enumerate(string):
	    if c == opened:
	        stack.append(i)
	    elif c == closed and stack:
	        start = stack.pop()
	        yield (len(stack), string[start + 1: i])

#********************
#get content in parenteses
def getParentheticContents(string, parenteses):
	return list(parentheticContents(string, parenteses))

#********************
#simplify KROME limits
def simplifyLimits(name):
	if name:
		name = name.replace(".LE.","<=").replace(".GE.",">=")
		name = name.replace(".LT.","<").replace(".GT.",">")
	return name

#********************
#limits to LaTeX format
def limits2latex(name):
	if name:
		name = name.replace("<=", " \leqslant ").replace(">="," \geqslant ")
		name = name.replace("K","\, \mathrm{K}")
	return name

#********************
#wrap exponential notation numbers in \num in latex string
def exp2latex(string):
        if(string):
                #replace 1e<x> by e<x> to get 10^x instead of 1x10^x
                string = re.sub("1e([+]*[0-9]{1,})", r"e\1",string)
                string = re.sub("([0-9]*\.*[0-9]*e[+-]*[0-9]{1,})", r"\\num{\1}",string)
        return string

#********************
#replaces fortran variable name with <replacement> in string
def replaceFortranVar(varname, replacement, string):
        string = re.sub("^"+varname+"$", replacement, string)
        string = re.sub("([^a-zA-Z])"+varname+"$", r'\g<1>'+replacement, string)
        string = re.sub("^"+varname+"([^a-zA-Z0-9_])", replacement+r'\g<1>', string)
        string = re.sub("([^a-zA-Z])"+varname+"([^a-zA-Z0-9_])", r'\g<1>'+replacement+r'\g<2>', string)
        return string

#********************
#adds line breaks to latex equation
def breakLatexEquation(string):
        opened = '{'
        closed = '}'
        depth = 0
	breakChars = ['+', '-']
        result = ""
        previousWasOpenParan = False
        numlines = 1
        for i, c in enumerate(string):
                if c == opened:
                        depth += 1
                elif c == closed:
			depth -= 1
                	#break on breakChars unless inside block or after opening paranthesis
                if c in breakChars and depth==0 and not previousWasOpenParan:
			c = " \\\\ \n &" + c
                        numlines += 1
                result += c
                if not c in ["(", " "]: previousWasOpenParan = False
                if c == "(": previousWasOpenParan = True

        return result, numlines

#********************
#replaces \left and \right by \bigL and \bigR in latex equation
# to allow parantheses to span multiple lines
def replaceLeftRightbyBigLR(rate):
        rate = rate.replace("\\left", "\\Bigl")
        rate = rate.replace("\\right", "\\Bigr")
        return rate

#********************
#raises curly brackets around latex expressions that belong to operators
def raiseBracketsOnOperators(rate):
        # Solution by scary recursive regular expression
        # ?2 means 'recursively add paranthesized group #2 here'
        # See e.g. http://www.rexegg.com/regex-recursion.html
        rate = regex.sub(r'(\operatorname\{[a-zA-Z]{1,}\})(\{(([^{}]|(?2))*)\})', r'\1\3', rate)
        return rate

def string_from(expr, string):
        for i in xrange(len(string)):
                if string[i:i+len(expr)] == expr:
                        splitAt = i + len(expr)
                        return (string[:splitAt], string[splitAt:])
        return string, ""

def next_parenthesis(string, left="{", right="}"):
        depth = 0
        acc = ""
        pos = 0
        
        for c in string:
                if c==left: depth += 1
                if c==right: depth -= 1
                acc += c
                if depth==0: break
                pos += 1
        
        return (acc, string[pos+1:])

def skip_spaces(string):
        acc = ""
        pos = 0
        for c in string:
                if c != " ": break
                pos += 1
        return (string[:pos],string[pos:])

#********************
#replaces fractions with divisor or dividend longer than maxlen by / sign
def replaceLongFracByDivide(rate,maxlen=100):
        acc=""
        while True:
                before, fromFrac = string_from("\\frac", rate)
                if len(fromFrac) <= 0:
                        return acc + rate
                spaces, rest = skip_spaces(fromFrac)
                numerator,rest = next_parenthesis(rest)
                denominator, rest = next_parenthesis(rest)
                
                rate = rest
                if max(len(numerator), len(denominator)) > maxlen:
                        acc += before[:-5] + spaces
                        acc += "\\left("+numerator[1:-1]+"\\right) / \\left("+denominator[1:-1] + "\\right)"
                else:
                        acc += before + spaces
                        acc += numerator + denominator

#********************
#resolve Krome variables in 'shortcuts', exept those in 'exceptions')
def replaceShortcuts(string, shortcuts, exceptions):
        for var in reversed(shortcuts):
                if var[0] in exceptions: continue
		if var[0] in string:
                        string = replaceFortranVar(var[0], var[1], string)
        return string

#********************
#get dictionary with shorcuts whose evaluation is deferred to end of latex table
def getDeferredShortcuts():
        from options import latexoptions
        return latexoptions.deferred_substitutions

#********************
#get table of symbols for latex table
def getSymbols():
        from options import latexoptions

        symbols = latexoptions.symbols.copy()

        # Add deferred symbols
        ds = getDeferredShortcuts()
        for key, value in ds.iteritems():
                symbols[key] = value

        return symbols

#********************
#get table of sympy symbols for latex table
def getSymbolTable():
        import sympy as sp

        symbols = getSymbols()

        # Convert to sympy symbols
        symboltable = {}
        for key in symbols:
                symboltable[key] = sp.Symbol(symbols[key])
        return symboltable


