from subprocess import call
import os,sys,hashlib,glob,platform

tests = ["collapse", "collapseZ", "dust", "map", "shock1D"]
tests += ["shock1Dphoto", "cloud","collapseUV","compact","lotkav","reverse"]
tests += ["shock1Dcool","slowmanifold"]

first = "" #start from this test (empty=first_test)
if(first.strip()==""): first = tests[0]

hashtab = [["collapse", "fort.22", "db231d69fbc3b04d1d4bb999aedd2c3d"],
	["collapseZ", "fort.55", "f8a9b722b8b7291bb9b641ac2e28a1e4"],
	["collapseZ", "fort.22", "885e29b3c08ca3692eefeb64ac0c0dc9"],
	['dust', 'fort.77', '5498df2824669fc98982f580d956433e'],
	['dust', 'fort.79', 'f1208a3f4e846a5a58aec050e8a8b3af'],
	['dust', 'fort.66', '800df2f862d3a315596c92401ad21b44'],
	['dust', 'fort.78', 'f7ac73aa11dbd45d79af963d6a64bf28'],
	['dust', 'fort.55', 'f8a9b722b8b7291bb9b641ac2e28a1e4'],
	['dust', 'fort.80', 'f1208a3f4e846a5a58aec050e8a8b3af'],
	['map', 'fort.66', '5bd3723abe293c599bd7b3af3f11dc1f'],
	['shock1D', 'fort.33', '8bcc5e0f6143f9458c7626f6ccacf98e'],
	['shock1D', 'fort.34', '928d1168f09cf27d74870f9f2177af8f'],
	['shock1D', 'fort.23', '3f1eb48618af4688564970d4bdbd1abe'],
	['shock1D', 'fort.24', '015062340b74694bebd119a0409d322e'],
	['shock1D', 'fort.22', '0443c7da947ac456dcc9b8e7b9a8572b'],
	["shock1Dphoto", "fort.33", "69bc834f5b16a299521241b810718c1e"],
	["shock1Dphoto", "fort.34", "7e2944d1e8589df91bf77027f54cf928"], 
	["shock1Dphoto", "fort.23", "ab2bad567db59a04a841e976e094e03e"],
	["shock1Dphoto", "fort.24", "c562df84de87addb753b5339e268b4ca"],
	["shock1Dphoto", "fort.22", "0443c7da947ac456dcc9b8e7b9a8572b"],
	["cloud", "fort.66", "ed6a8a52bec401d3167c8640e2171a72"],
	["collapseUV", "fort.22", "0dfcdedf30849458dc8c6df9b02de962"],
	['compact', 'fort.33', '62a8dd458361aab8fdcf8067d530ee9c'],
	['compact', 'fort.34', '928d1168f09cf27d74870f9f2177af8f'],
	['compact', 'fort.23', '3f1eb48618af4688564970d4bdbd1abe'],
	['compact', 'fort.24', '015062340b74694bebd119a0409d322e'],
	['compact', 'fort.22', '0443c7da947ac456dcc9b8e7b9a8572b'],
	["lotkav", "fort.66", "a560834860a7045c71842831d6765d8f"],
	["reverse", "fort.66", "0696206b1578294e14485e167f36d757"],
	['shock1Dcool', 'fort.33', '4b3faafd1b37f14b4a4cbe7a2a55d01a'],
	['shock1Dcool', 'fort.34', 'a2ddf647cf76a51b85c7132e0c010827'],
	['shock1Dcool', 'fort.23', 'bcb9fffb39f991d437bd508335be81b4'],
	['shock1Dcool', 'fort.24', '8829dc9c235e4009c7b11d2d8a0d9cd0'],
	['shock1Dcool', 'fort.22', '0443c7da947ac456dcc9b8e7b9a8572b'],
	["slowmanifold", "fort.77", "d5854cf7627dc3f29c0213256d933a20"],
	["slowmanifold", "fort.66", "a1fdbcab746107c4c834cc84a1f7c06b"]]

print "************************************************"
print "WARNING: this script will ERASE all the contents"
print " in the ./build folder!"
print "************************************************"
a = raw_input("Any key to continue q to quit... ")
if(a=="q"): print sys.exit()

print
print "************************************************"
print "WARNING: this tests run with -O0 option and all"
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
	#hash and store MD5 for fort files
	hashall = []
	for fort in glob.glob("fort.*"):
		md5 = hashlib.md5(open(fort).read()).hexdigest()
		hashall.append([test,fort,md5])
		fout.write(test+" "+fort+" "+md5+"\n")
	for hashblock in hashall:
		if(hashblock not in hashtab):
			print "ERROR with "+(",".join(hashblock))
			for x in hashall:
				print x
			sys.exit()

	#gnuplot if you want graphical result
	#call(["gnuplot"])
	#clear directory
	for ff in glob.glob("*"):
		os.remove(ff)
	#back to krome main directory
	os.chdir("..")
	print "DONE!"


	
	

