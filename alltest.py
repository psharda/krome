from subprocess import call
import os,sys,hashlib,glob,platform

tests = ["collapse", "collapseZ", "dust", "map", "shock1D"]
tests += ["shock1Dphoto", "cloud","collapseUV","compact","lotkav","reverse"]
tests += ["shock1Dcool","slowmanifold"]

first = "collapse" #start from this test

hashtab = [["collapse", "fort.22", "0a53c195e1343ba0248edce4a8611455"],
	["collapseZ", "fort.55", "f8a9b722b8b7291bb9b641ac2e28a1e4"],
	["collapseZ", "fort.22", "6df2c545d3dd1384e0a9c7514d240c80"],
	["dust", "fort.77", "9965fc66f87c37c696cfefd534ba7c8b"],
	["dust", "fort.79", "f1208a3f4e846a5a58aec050e8a8b3af"],
	["dust", "fort.66", "1a7e27a5b115fa7a2b14a946e815983d"],
	["dust", "fort.78", "b2a82cb1926c309cf876899320b0d1d6"],
	["dust", "fort.80", "f1208a3f4e846a5a58aec050e8a8b3af"],
	["map", "fort.66", "5bd3723abe293c599bd7b3af3f11dc1f"],
	["shock1D", "fort.33", "225ec1cda13d44ef67d932cdea96ee33"],
	["shock1D", "fort.34", "73d25aeb344de4b146512c94b9349673"],
	["shock1D", "fort.23", "3f1eb48618af4688564970d4bdbd1abe"],
	["shock1D", "fort.24", "015062340b74694bebd119a0409d322e"],
	["shock1D", "fort.22", "0443c7da947ac456dcc9b8e7b9a8572b"],
	["shock1Dphoto", "fort.33", "69bc834f5b16a299521241b810718c1e"],
	["shock1Dphoto", "fort.34", "7e2944d1e8589df91bf77027f54cf928"],
	["shock1Dphoto", "fort.23", "ab2bad567db59a04a841e976e094e03e"],
	["shock1Dphoto", "fort.24", "c562df84de87addb753b5339e268b4ca"],
	["shock1Dphoto", "fort.22", "0443c7da947ac456dcc9b8e7b9a8572b"],
	["cloud", "fort.66", "ed6a8a52bec401d3167c8640e2171a72"],
	["collapseUV", "fort.22", "0dfcdedf30849458dc8c6df9b02de962"],
	["compact", "fort.33", "9291a9a1ab23d1a7f56bf865d1c0459b"],
	["compact", "fort.34", "73d25aeb344de4b146512c94b9349673"],
	["compact", "fort.23", "3f1eb48618af4688564970d4bdbd1abe"],
	["compact", "fort.24", "015062340b74694bebd119a0409d322e"],
	["compact", "fort.22", "0443c7da947ac456dcc9b8e7b9a8572b"],
	["lotkav", "fort.66", "a560834860a7045c71842831d6765d8f"],
	["reverse", "fort.66", "0696206b1578294e14485e167f36d757"],
	["shock1Dcool", "fort.33", "7e966d8124ad14da937731e9e688b322"],
	["shock1Dcool", "fort.34", "1a0df9be88904e626d1b5074c2bbd0aa"],
	["shock1Dcool", "fort.23", "2c62760bd667e6810baf1715137469f6"],
	["shock1Dcool", "fort.24", "4cb4133fa03e81ca42c4a4d7661a0ead"],
	["shock1Dcool", "fort.22", "0443c7da947ac456dcc9b8e7b9a8572b"],
	["slowmanifold", "fort.77", "d5854cf7627dc3f29c0213256d933a20"],
	["slowmanifold", "fort.66", "a1fdbcab746107c4c834cc84a1f7c06b"]]

print "************************************************"
print "WARNING: this script will ERASE all the contents"
print " in the ./build folder!"
print "************************************************"
a = raw_input("Any key to ignore q to quit... ")
if(a=="q"): print sys.exit()

#open output file for MD5
fout = open("outtest","w")
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
			call([notifier,"--notification","--text", "\""+test+" done!\"", "--timeout","2"])
	#hash and store MD5 for fort files
	for fort in glob.glob("fort.*"):
		md5 = hashlib.md5(open(fort).read()).hexdigest()
		hashblock = [test,fort,md5]
		if(hashblock not in hashtab):
			print "ERROR with "+(",".join(hashblock))
			sys.exit()
		fout.write(test+" "+fort+" "+md5+"\n")
	#gnuplot if you want graphical result
	#call(["gnuplot"])
	#clear directory
	for ff in glob.glob("*"):
		os.remove(ff)
	#back to krome main directory
	os.chdir("..")
	print "DONE!"


	
	

