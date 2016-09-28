from subprocess import call
from shutil import copyfile
import ftplib,os,glob,hashlib

makeOption = "debug" #option makefile
testpath = "tests/" #where the tests are located
prj_name = "alltest" #where to test
plotFolder = "plot_alltests" #plots folder
pullKp = False

#move to kp folder
os.chdir("../kp/")
if(pullKp): call(["git", "pull","origin"])

#go back to krome folder
os.chdir("../krome/")

#import the list of tests from testpath
tests = [x[0].replace(testpath,"") for x in os.walk(testpath) if x[0]!=testpath]

#create plot directory if none
if(not(os.path.exists(plotFolder))):
	os.makedirs(plotFolder)

#clear plot directory
for ff in glob.glob(plotFolder+"/*"):
	os.remove(ff)


hashall = []

#loop on tests
for test in tests[0:3]:
	print
	print "#########################################################################################################"
	print "                                          test: "+test
	print "#########################################################################################################"
	print

	#clear directory
	for ff in glob.glob("build_"+prj_name+"/*"):
		os.remove(ff)

	#call krome
	callarg = ["./krome", "-test="+test,"-skipDevTest", "-unsafe","-project="+prj_name]
	print callarg
	call(callarg)

	#skip development test
	if(not(os.path.isfile("build_"+prj_name+"/Makefile"))): continue

	#move to build directory
	os.chdir("build_"+prj_name+"/")

	#make clean
	call(["make","clean"])

	#compile full debug
	call(["make",makeOption])

	#run executable
	call(["./test"])

	#copy kp class to current build folder
	copyfile("../../kp/kp.py", "kp.py")

	#make clean
	call(["python","plot.py"])

	#copy plot to plot folder
	copyfile("plot.png", "../"+plotFolder+"/"+test+".png")

	#hash and store MD5 for fort files
	for fort in glob.glob("fort.*"):
		md5 = hashlib.md5(open(fort).read()).hexdigest()
		hashall.append([makeOption,test,fort,md5])

	#move to build directory
	os.chdir("..")

fout = open("hash.result","w")
#loop on md5 hashes
for xhash in hashall:
	print xhash
	fout.write((" ".join(xhash))+"\n")
fout.close()



