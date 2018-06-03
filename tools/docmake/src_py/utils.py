import os,datetime,math

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
def char2int(arg):
	if isNumber(arg):
		return int(float(arg))
	else:
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
	shortcut = []
	fileName = "temperatureShortcuts.dat"
        absPath = os.path.join(os.path.dirname(__file__), "..", fileName)
        absPath = os.path.abspath(absPath)
	fh = open(absPath,"rb")
	for row in fh:
		srow = row.strip()
		if(srow==""): continue
		if(srow.startswith("#")): continue
		(variable,expression) = [x.strip() for x in srow.split("!")[0].split("=")]
		shortcut.append((variable,expression))
	return shortcut

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
