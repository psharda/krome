import sys
#this python script returns the list and the information on the functions
# in the krome_user.f90 file
# the details are taken from the block of comments before the function when present
# Grassi+KROME Team 11Jul2014

#NO NEED TO MODIFY THIS SCRIPT
listNamesOnly = False
argv = sys.argv
if(len(argv)>1):
	if("-h" in argv):
		print "usage: python "+argv[0]+" [OPTIONS]"
		print "\t-n list names only"
		print "\t-h show this help and exit"
		sys.exit()
	listNamesOnly = ("-n" in argv)
		

fh = open("krome_user.f90","rb")

inpartF = ["function","(",")"]
inpartS = ["subroutine","(",")"]
infun = issub = False
storecom = ""
flist = []
for row in fh:
	srow = row.strip()
	if(srow==""): continue
	if(srow=="contains"): storecom = ""
	isfun = True
	for part in inpartF:
		if(not(part in srow)): 
			isfun = False
			break
	issub = True
	for part in inpartS:
		if(not(part in srow)): 
			issub = False
			break
	if(isfun or issub):
		if(storecom==""): storecom = "[no comments available]"
		ffname = srow
		flist.append([ffname,storecom.strip()])
		continue

	if(("end function" in srow) or ("end subroutine" in srow)):
		storecom = ""

	if(srow[0]=="!" and not("*******" in srow)): storecom += "  "+srow[1:].strip()+"\n"

flist = sorted(flist,key=lambda x:x[0])
icount = 0
for x in flist:
	print str(icount+1)+") "+x[0]
	print "  "+x[1]+"\n"
	icount += 1

