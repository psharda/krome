#this script runs all tests and/or produces or checks md5
from subprocess import call
import os,sys,hashlib,glob,platform
tests = ["atmosphere", "auto","chianti","collapseCO","collapseUV_Xrays","collapseZ_induced"]
tests += ["collapse", "collapseZ", "dust", "map", "shock1D","compact","earlyUniverse","lamda"]
tests += ["cloud","collapseUV","compact","lotkav","reverse","wrapC"]
tests += ["shock1Dcool","slowmanifold"]

first = "lamda" #start from this test (empty string start from first test)
if(first.strip()==""): first = tests[0]

mode = "eyeball" #"hash":produce hashfile, "eyeball":hashfile+call gnuplot to plot ,"": check hash

#read hastable if needed
if(mode=="hash"):
	hashtab = []
	fh = open("outtest.log","rb")
	for row in fh:
		srow = row.strip()
		arow = srow.split(" ")
		hashtab.append(arow)

print "************************************************"
print "WARNING: this script will ERASE all the contents"
print " in the ./build folder!"
print "************************************************"
a = raw_input("Any key to continue q to quit... ")
if(a=="q"): print sys.exit()

print
print "************************************************"
print "WARNING: these tests run with -O0 option and all"
print " the check flags enabled, so they are very slow!"
print "************************************************"
a = raw_input("Any key to continue q to quit... ")
if(a=="q"): print sys.exit()


os.chdir("build/")
#clear directory
for ff in glob.glob("*"):
	os.remove(ff)
os.chdir("..")

#open output file for MD5
if(mode!=""): 
	fout = open("outtest.log","w")
run = False #run flag
for test in tests:
	if(test==first): run = True #run the first test
	if(not(run)): continue #skip if test is before first
	print "test "+test
	#call krome
	call(["./krome", "-test="+test, "-pedantic", "-unsafe", "-sh"])
	#move to build directory
	os.chdir("build/")
	#make clean
	call(["make","clean"])
	#compile
	call(["make"])
	#run executable
	call(["./krome"])
	#run zenity notification when exectutable ends (LINUX USERS ONLY)
	if("linux" in platform.system().lower()):
		notifier = "zenity"
		fpath = "/usr/bin/"+notifier
		#check if notifier zenity exists
		if(os.path.isfile(fpath) and os.access(fpath, os.X_OK)):
			call(["killall","notification-daemon"])
			call([notifier,"--notification","--text", "\""+test+" done!\"", "--timeout","2"])

	hashall = []
	#hash and store MD5 for fort files
	for fort in glob.glob("fort.*"):
		md5 = hashlib.md5(open(fort).read()).hexdigest()
		hashall.append([test,fort,md5])
		if(mode!=""): fout.write(test+" "+fort+" "+md5+"\n")
	#control the hash found
	if(mode==""):
		for hashblock in hashall:
			if(hashblock not in hashtab):
				print "ERROR with "+(",".join(hashblock))
				for x in hashall:
					print x
				sys.exit()

	#gnuplot if you want graphical result
	if(mode=="eyeball"): call(["gnuplot"])
	#clear directory
	for ff in glob.glob("*"):
		os.remove(ff)
	#back to krome main directory
	os.chdir("..")
	print "DONE!"


	
	

