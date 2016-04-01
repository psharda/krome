#!/usr/bin/python
#!python

#This script print the missing reverse reactions with an arbitrary small value
#written by Tommaso Grassi and the KROME team (Apr,2016)

import sys,os
from subprocess import call
from random import random as rand
from math import sin,cos

########################################
# PLEASE DO NOT MODIFY IF YOU ARE NOT 
# SURE OF WHAT YOU ARE DOING!
########################################

#check command-line arguments and store them
if(len(sys.argv)<1):
	print ("Usage: %s INPUT" % sys.argv[0])
	sys.exit()

if(not os.path.isfile(sys.argv[1].strip())):
    sys.exit("ERROR: input file %s was not found!" % sys.argv[1])


smallValue = 0e0
fname = sys.argv[1].strip()

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

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

cfmt = form() #object class format
fh = open(fname,"rb") #open file to read
foundrea = [] #list of reaction found (for output on screen only)
foundrev = [] #list of reaction found (for output on screen only)


print "**************************"
print "        SMALL REVERSE!"
print "**************************"
print "Reading from "+fname
#default format
fmt = "@format:idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"
cfmt.proc(fmt) #process format
skip = False
ifound = 0
for row in fh:
	#skip comments and blanks
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
		continue
	if(skip): continue

	#process or skip format statement
	if("@format:" in srow):
		cfmt.proc(srow)
		continue
	if(srow.startswith("@")): continue
	if(srow.startswith("if")): continue

	#split line
	arow = srow.split(",")

	#find reactants
	rr = []
	for ir in cfmt.reactants:
		if(arow[ir].strip()!=""): rr.append(arow[ir])

	#find products
	pp = []
	for ip in cfmt.products:
		if(arow[ip].strip()!=""): pp.append(arow[ip])
	foundrea.append([sorted(rr),sorted(pp),"fwd"])
	#add reverse
	foundrea.append([sorted(pp),sorted(rr),"rev"])

#check if reverse is present already as forward
block = []
oldRea = [None]*3
revCount = 500
fmtOld = ""
#loop on found reactions
for rea in sorted(foundrea):
	if(oldRea[:2]==rea[:2]):
		block.append(rea)
	else:
		if(block!=[]):
			#check if forward if present
			fwdFound = False
			for b in block:
				if(b[2]=="fwd"): fwdFound = True
			#if fwd is not present print reverse
			if(not(fwdFound)):
				brea = block[0]
				fmt = "@format:idx,Tmin,Tmax,"+("R,"*len(brea[0]))+("P,"*len(brea[1]))+"rate"
				if(fmt!=fmtOld): print fmt
				print str(revCount)+",NONE,NONE,"+(",".join(brea[0]))+","+(",".join(brea[1]))+","+str(smallValue)
				revCount += 1
				fmtOld = fmt

		block = [rea]
	oldRea = rea


print "DONE, bye!"



