#this python script returns the list and the information on the functions
# in the krome_user.f90 file
# the details are taken from the block of comments before the function when present
# Grassi+KROME Team 11Jul2014


#NO NEED TO MODIFY THIS SCRIPT
fh = open("krome_user.f90","rb")

inpartF = ["function","(",")"]
inpartS = ["subroutine","(",")"]
infun = issub = False
storecom = ""
icount = 0
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
		if(storecom==""): storecom = "  [no comments available]\n"
		icount += 1
		print str(icount)+") "+srow
		print storecom
		continue

	if(("end function" in srow) or ("end subroutine" in srow)):
		storecom = ""

	if(srow[0]=="!" and not("*******" in srow)): storecom += "  "+srow[1:].strip()+"\n"

