#this script runs all tests and/or produces or checks md5
#it is employed to run automatic tests
#no need to run or modify by the user
from subprocess import call
from subprocess import Popen, PIPE
import os,sys,hashlib,glob,platform,shutil,time
import ftplib,urllib2,sys,cookielib

argv = sys.argv

compiler="ifort" #ifort or gfortran

if("-compiler" in argv):
	compiler = argv[argv.index("-compiler")+1]

testpath = "tests/" #where the tests are located
prj_name = "alltest" #where the tests

#import the list of tests from testpath
tests = [x[0].replace(testpath,"") for x in os.walk(testpath) if x[0]!=testpath]
#tests = ["hello"]

#start from this test (empty string start from first test)
first = ""
if(first.strip()==""): first = tests[0]


mode = "check" #"hash":produce hashfile, "eyeball":hashfile+call gnuplot to plot ,"check": check hash

################

#check if a new test is needed
if(mode=="check"):
	#synchronize
	call(["git", "pull", "origin"])

	changesetFOLDER = "" #changeset of the current folder
	proc = Popen(['git','show'],stdout=PIPE) #use git show to retrive local info
	#loop on the output
	for line in iter(proc.stdout.readline,''):
		lstrip = line.rstrip()
		#grep the commit line
		if("commit" in lstrip):
			changesetFOLDER = lstrip.replace("commit","").strip()[:7]
			break
	#check if the changeset is retrieved
	if(changesetFOLDER==""): sys.exit("ERROR: check mode enabled and git show command does not work properly")

	#retireve the changeset on the SERVER
	changesetSERVER = ""
	#content = urllib2.urlopen('http://kromepackage.org/test/outcheck.log')
	cj = cookielib.CookieJar()
	opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
	request = urllib2.Request('http://kromepackage.org/test/outcheck.log')
	content = opener.open(request)
	for line in iter(content.readlines()):
		if("changeset:" in line): changesetSERVER = line.replace("changeset:","").strip()
	#check if the server changeset is retrieved
	if(changesetSERVER==""): sys.exit("ERROR: check mode enabled and url not retrieved")

	#check if the server and the folder have the same changeset. If so no need to check.
	#if(changesetSERVER==changesetFOLDER):
	#	sys.exit("SERVER and local FOLDER has the same changeset. No need to check.")



#read hastable if needed
if(mode=="check"):
	hashtab = []
	changeset = "unknown"
	if(not(os.path.isfile("outtest.md5"))): sys.exit("ERROR: file outtest.md5 not present!")
	fh = open("outtest.md5","rb")
	for row in fh:
		srow = row.strip()
		#load reference changeset name
		if("changeset:" in srow):
			changesetREF = srow.replace("changeset:","")
			continue
		arow = srow.split(" ")
		hashtab.append(arow)

#print "************************************************"
#print "WARNING: this script will ERASE all the contents"
#print " in the ./build folder!"
#print "************************************************"
#a = raw_input("Any key to continue q to quit... ")
#if(a=="q"): print sys.exit()

#print
#print "************************************************"
#print "WARNING: these tests run with -O0 option and all"
#print " the check flags enabled, so they are very slow!"
#print "************************************************"
#a = raw_input("Any key to continue q to quit... ")
#if(a=="q"): print sys.exit()

if(not(os.path.exists("build_"+prj_name+"/"))): os.makedirs("build_"+prj_name+"/")

os.chdir("build_"+prj_name+"/")
#clear directory
for ff in glob.glob("*"):
	os.remove(ff)
os.chdir("..")

#get changeset
masterfile = ".git/refs/heads/master" #name of the git master file
changeset = "unknown"
#if git master file exists grep the changeset
if(os.path.isfile(masterfile)):
	changeset = open(masterfile,"rb").read()
	changeset = changeset[:7]

#open output file for MD5 or results
if(mode!="check"):
	fout = open("outtest.log","w")
else:
	fout = open("outcheck.log","w")
	fout.write("changesetREF: "+changesetREF+"\n")

#write the changeset to file
fout.write("changeset: "+changeset+"\n")

run = False #run flag
for test in tests:
	if(test==first): run = True #run the first test
	if((compiler=="gfortran") and (test=="wrapC")): continue
	if(not(run)): continue #skip if test is before first
	print
	print "#########################################################################################################"
	print "                                          test: "+test
	print "#########################################################################################################"
	print

	#clear directory
	for ff in glob.glob("build_"+prj_name+"/*"):
		os.remove(ff)

	#call krome
	callarg = ["./krome", "-test="+test,"-skipDevTest", "-pedantic", "-unsafe", "-sh","-compiler="+compiler,"-project="+prj_name]
	print callarg
	call(callarg)

	#skip development test
	if(not(os.path.isfile("build_"+prj_name+"/Makefile"))): continue

	#move to build directory
	os.chdir("build_"+prj_name+"/")

	#make clean
	call(["make","clean"])

	#compile
	call(["make"])

	#run executable
	call(["./krome"])

	#run zenity notification when exectutable ends (LINUX USERS ONLY)
	#if("linux" in platform.system().lower()):
	#	notifier = "zenity"
	#	fpath = "/usr/bin/"+notifier
	#	#check if notifier zenity exists
	#	if(os.path.isfile(fpath) and os.access(fpath, os.X_OK)):
	#		call(["killall","notification-daemon"])
	#		call([notifier,"--notification","--text", "\""+test+" done!\"", "--timeout","2"])

	hashall = []
	#hash and store MD5 for fort files
	for fort in glob.glob("fort.*"):
		md5 = hashlib.md5(open(fort).read()).hexdigest()
		hashall.append([test,fort,md5])
		print test,fort,md5
		if(mode!="check"): fout.write(test+" "+fort+" "+md5+"\n")

	#control the hash found if needed
	if(mode=="check"):
		print "checking..."
		testOK = True
		for hashblock in hashall:
			if(hashblock not in hashtab):
				print "ERROR with "+(",".join(hashblock))
				testOK = False
				for x in hashall:
					print x
				#sys.exit()
		print "Is OK?",testOK
		fout.write(test+" "+str(testOK)+" "+str(time.time())+" regular\n")


	#call gnuplot if you want graphical result
	if(mode=="eyeball"): call(["gnuplot"])

	#clear directory
	for ff in glob.glob("*"):
		os.remove(ff)

	#back to krome main directory
	os.chdir("..")
	print "DONE!"

#add skipped tests as failed
#if(mode=="check"):
#	for test in skiptests:
#		fout.write(test+" "+str(False)+" "+str(time.time())+" skipped\n")
fout.close()

#copy the results to kromepackage.org using FTP
import traceback
if(mode=="check"):
	filename = "outcheck.log"
	if(not(os.path.isfile(filename))): sys.exit("ERROR: "+filename+" not present. Nothing to copy.")
	if(not(os.path.isfile("../ftplogin.dat"))): sys.exit("ERROR: ftplogin.dat not present. Can't connect.")
	print "copying "+filename+" using FTP..."
	usr,psw = [x.strip() for x in open("../ftplogin.dat","rb").read().split()]
	ftp = ftplib.FTP("kromepackage.org",usr, psw)
	try:
		ftp.cwd("/test")
		ftp.storbinary("STOR "+filename, file(filename,"rb"))
	except:
		traceback.print_exc()
	ftp.quit()
	print "DONE!"

