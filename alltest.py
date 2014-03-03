from subprocess import call
import os,sys,hashlib,glob

tests = ["collapse", "collapseZ", "dust", "map", "shock1D"]
tests += ["shock1Dphoto", "cloud","collapseUV","compact","lotkav","reverse"]
tests += ["shock1Dcool","slowmanifold"]

first = "collapse" #start from this test

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
	#run zenity notification (LINUX USERS ONLY)
	call(["zenity","--notification","--text", "\""+test+" done!\"", "--timeout","2"])
	#hash and store MD5 for fort files
	for fort in glob.glob("fort.*"):
		fout.write(test+" "+fort+" "+hashlib.md5(open(fort).read()).hexdigest()+"\n")
	#gnuplot if you want graphical result
	call(["gnuplot"])
	#clear directory
	for ff in glob.glob("*"):
		os.remove(ff)
	#back to krome main directory
	os.chdir("..")
	print "DONE!"


	
	

