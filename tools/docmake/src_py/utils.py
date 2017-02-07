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

#***********************
def getFooterInfo():
	#name of the git master file
	masterfile = "../../.git/refs/heads/master"
	changeset = ("x"*7) #default unknown changeset
	#if git master file exists grep the changeset
	if(os.path.isfile(masterfile)):
		changeset = open(masterfile,"rb").read()

	datenow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	bitbucketLink = "https://bitbucket.org/tgrassi/krome/commits/"+changeset
	bitbucketLanding = "<a href=\"https://bitbucket.org/tgrassi/krome\" target=\"_blank\">KROME</a>"
	hrefWiki = "<a href=\"https://bitbucket.org/tgrassi/krome/wiki/docmaker\" target=\"_blank\">docmake</a>"
	return "documentation generated with "+hrefWiki+" ("+bitbucketLanding+") - changeset: <a href=\""\
		+ bitbucketLink + "\" target=\"_blank\">" + changeset[:7] + "</a> - " + datenow

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
	if(arg==0): return "0"
	if(arg>=maxLim):
		return htmlExp(arg,digits=digits)
	else:
		return str(arg)

#*********************
def htmlExp(arg,digits=2):
	if(arg==0): return "0"
	xp = int(math.log10(abs(arg)))
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

