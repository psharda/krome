#!/usr/bin/python
#!python

import sys,os
#THIS UTILITY CHECKS A RATE FILE
#1) CHECK FOR DUPLICATES
#2) CHECK FORMAT
#3) CHECK BRACKETS
#4) RENUMBER RATES
#written by Tommaso Grassi and the KROME team (Mar21,2014)

#check command-line arguments and store them
if(len(sys.argv)<3):
    sys.exit("Usage: %s INPUT OUTPUT" % sys.argv[0])

if(not os.path.isfile(sys.argv[1].strip())):
    sys.exit("ERROR: input file %s was not found!" % sys.argv[1])

fname = sys.argv[1]
outfile = sys.argv[2]
if(fname==outfile):
    sys.exit("ERROR: input and output are the same file!")

##############
class form():
	reactants = []
	products = []
	idx = 0
	tmin = 0
	tmax = 0
	rate = 0
	fmt = ""
	nfmt = 0
	def proc(self,myfmt):
		self.reactants = []
		self.products = []
		self.idx = 0
		self.tmin = 0
		self.tmax = 0
		self.rate = 0
		self.fmt = myfmt
		self.nfmt = len(myfmt.split(","))
		afmt = myfmt.replace("@format:","").lower().split(",")
		while("r" in afmt):
			ii = afmt.index("r")
			self.reactants.append(ii)
			afmt[ii] = ""
		while("p" in afmt):
			ii = afmt.index("p")
			self.products.append(ii)
			afmt[ii] = ""
		if("idx" in afmt):
			self.idx = afmt.index("idx")
		if("tmin" in afmt):
			self.tmin = afmt.index("tmin")
		if("tmax" in afmt):
			self.tmax = afmt.index("tmax")
		if("rate" in afmt):
			self.rate = afmt.index("rate")
	def show(self):
			print self.reactants
			print self.products
			print self.idx
			print self.tmin
			print self.tmax
			print self.rate
			print self.fmt

okall = True
uniq = []
cfmt = form()
fh = open(fname,"rb")
fout = open(outfile,"w")
icount = 0
print "**************************"
print "    REACTION CHECKER"
print "**************************"
print "Reading from "+fname
fmt = "@format:idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"
cfmt.proc(fmt)
skip = False
wasComment = True
for row in fh:
	skipThis = False
	srow = row.strip()
	if(srow==""): continue
	if(srow[0]=="#"): skipThis = True
	if(srow[0:2]=="//"): skipThis = True
	if("*/" in srow):
		skip = False
		skipThis = True
	if(srow[0:2]=="/*"): 
		skip = True
		skipThis = True
	if(skipThis):
		if(not(wasComment)): fout.write("\n")
		wasComment = True
		fout.write(srow+"\n")
		continue
	if(skip): continue
	#print srow

	if(srow[0]=="@"): fout.write(srow+"\n")

	if("@format:" in srow):
		cfmt.proc(srow)
		continue
	if(srow[0]=="@"): continue

	wasComment = False
	arow = srow.split(",",cfmt.nfmt-1)
	if(len(arow)<cfmt.nfmt):
		print "***********************"
		print "ERROR: wrong format for"
		print srow
		print "Expected "+str(cfmt.nfmt)+" parts, found "+str(len(arow))
		okall = False
		print cfmt.fmt
		sys.exit()

	if(len(arow)>cfmt.nfmt):
		print "***********************"
		print "WARNING: possible wrong format for"
		print srow
		print "Expected "+str(cfmt.nfmt)+" parts, found "+str(len(arow))
		okall = False
		print cfmt.fmt

	if(arow[cfmt.rate].count("(")!=arow[cfmt.rate].count(")")):
		print "***********************"
		print "ERROR: unbalanced brackets!"
		print srow
		sys.exit()

	rr = []
	for ir in cfmt.reactants:
		if(arow[ir].strip()!=""): rr.append(arow[ir])
	pp = []
	for ip in cfmt.products:
		if(arow[ip].strip()!=""): pp.append(arow[ip])
	tt = ["",""]
	if(cfmt.tmin!=0): tt[0] = arow[cfmt.tmin].lower().replace("none","")
	if(cfmt.tmax!=0): tt[1] = arow[cfmt.tmax].lower().replace("none","")
	rr = sorted(rr)
	pp = sorted(pp)
	for uu in uniq:
		if([rr,pp,tt] == uu[0]):
			print "***********************"
			print "WARNING: reaction already found"
			print rr,pp,tt
			print srow
			print uu[1]
			okall = False
			#cfmt.show()

	uniq.append([[rr,pp,tt],srow])
	icount += 1
	arow[cfmt.idx] = str(icount)
	fout.write((",".join(arow))+"\n")
print str(icount)+" reactions found"
if(okall): print "Everything seems ok!"
print "Output written in "+outfile
fout.close()
print "DONE, bye!"

