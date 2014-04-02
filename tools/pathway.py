#!/usr/bin/python
#!python

#THIS CREATES A GRAPH WITH GRAPHVIZ
# FOR ALL THE SPECIES THAT CONTAINS A GIVEN ATOM
#written by Tommaso Grassi and the KROME team (Mar25,2014)

import sys,os
from subprocess import call
from random import random as rand
from math import sin,cos

########################################
# PLEASE DO NOT MODIFY IF YOU ARE NOT 
# SURE OF WHAT YOU ARE DOING!
########################################
tmpfile = str(1e8+rand()*1e8)+"tmp.dot"
EMAX = 2 #maximum number of object on label edge. if greater select by degree
DEGMIN = 0 #limit on degree percentage nodes
DEGMIN_LABEL = 0 #limit on degree percentage edge
prog = "dot" #graphviz tool employed



#check command-line arguments and store them
if(len(sys.argv)<4):
	print ("Usage: %s INPUT OUTPUT BASE" % sys.argv[0])
   	print " where BASE is the atom shared by the species (case sensitive!)"
	sys.exit()

if(not os.path.isfile(sys.argv[1].strip())):
    sys.exit("ERROR: input file %s was not found!" % sys.argv[1])

fname = sys.argv[1].strip()
outfile = sys.argv[2].strip()
base = sys.argv[3].strip()
if(fname==outfile):
    sys.exit("ERROR: input and output are the same file!")
if(not(".eps" in outfile)):
    sys.exit("ERROR: output file should end with .eps!")

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

adic = ['H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL','SI','P','S','CL','AR',\
		'K','CA','SC','TI','V','CR','MN','FE','NI','CU','ZN','GA','GE','AS','SE','BR','KR',\
		'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE',\
		'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',\
		'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',\
		'AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U']
adic = sorted([x.lower() for x in adic],key=lambda x:len(x),reverse=True)

def basein(mydic,argin,hayin):
	hay = hayin.lower()
	arg = argin.lower()
	for el in mydic:
		if(len(el)>len(arg)): hay = hay.replace(el,"")
	return (arg in hay)

cfmt = form() #object class format
fh = open(fname,"rb") #open file to read
fout = open(tmpfile,"w") #open output file
foundedge = dict() #dictionary key="x -> y", value=[reaction partners]
foundrea = [] #list of reaction found (for output on screen only)

#start writin digraph file
fout.write("//DOT file created with PATHWAY tool from KROME\n")
fout.write("// see kromepackage.org\n")
fout.write("digraph pathway{\n")
fout.write("node [shape=circle];\n")
#fout.write("nodesep=0.6;\n")
fout.write("splines=true;\n")
fout.write("overlap=scalexy;\n")
#fout.write("sep=\"+25,25\";\n")

fout.write("fontsize=11;\n")

degree = dict()

print "**************************"
print "        PATHWAYS!"
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
	if(srow[0]=="@"): continue

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
	#check for reaction with more than 3.body
	if(len(rr)>3): 
		print "WARNING: "+str(len(rr))+"-body skipped!"
		continue

	#store degree
	crhv = []
	if(len(rr)==1): crhv = ["CR/hv"] #empty reactants means CR or photon
	for x in rr+pp+crhv:
		if(x in degree): 
			degree[x] += 1
		else:
			degree[x] = 1

	#search for the atom in reactants and products
	rrpp = "".join(rr+pp)
	if(basein(adic,base,rrpp)):
		irr = 0
		#loop on reactants
		for x in rr:
			#if a reactants has the atom
			if(basein(adic,base,x)):
				node_start = x #store reactants
				other_reacts = []
				#store all other reactants
				for i in range(len(rr)):
					if(i!=irr): other_reacts.append(rr[i])

				node_ends = []
				#store products with the atom
				for p in pp:
					if(basein(adic,base,p)):
						node_ends.append(p)
				#empty reactants are CR or photons
				if(len(other_reacts)==0): other_reacts = ["CR/hv"]
				#loop on products
				for nodee in node_ends:
					ss = "\""+x+"\" -> \""+nodee+"\"" #edge as key
					#add reactants
					if(ss in foundedge):
						foundedge[ss] += other_reacts
					else:
						foundedge[ss] = other_reacts

					#check for reactions to print on screen
					if([rr,pp] in foundrea): continue
					ifound += 1
					print str(ifound)+". "+(" + ".join(rr))+" -> "+(" + ".join(pp))
					foundrea.append([rr,pp])

			irr += 1 #count reactants

if(ifound==0):
	sys.exit("ERORR: no reactions found with atom "+base+"!")

#dictionary->list (easier to sort)
sumdeg = sum(degree.values())
deg = []
degnorm = dict()
for k,v in degree.iteritems():
	deg.append([k,v*100./sumdeg])
	degnorm[k] = int(v*100./sumdeg)
deg = sorted(deg, key=lambda x:x[1],reverse=True)
degree = degnorm
#print deg

#prepare dot file
species = []
gotwarn = False
for k,v in foundedge.iteritems():
	uniq = []
	#unique labels on edge
	for x in v:
		if(x in uniq): continue
		uniq.append(x)
	node_start = k.split("->")[0].strip().replace("\"","")
	node_end = k.split("->")[1].strip().replace("\"","")
	species.append(node_start)
	species.append(node_end)	
	if(degree[node_start]<DEGMIN or degree[node_end]<DEGMIN): continue
	#print "*********"
	#print k
	#print uniq
	maxdeg = 0
	#eliminate edge with partners with max degree < DEGMIN_LABEL
	for x in uniq:
		maxdeg = max(degree[x],maxdeg)
	if(maxdeg<DEGMIN_LABEL): continue

	vals = []
	#if more than EMAX objects on the edge label keeps only the EMAX with highes degree
	if(len(uniq)>EMAX):
		tmpu = (", ".join(uniq))
		for x in uniq:
			vals.append([x,degree[x]])
		vals = sorted(vals, key=lambda x:x[1],reverse=True)
		uniq = [x[0] for x in vals[0:EMAX]]
		print "WARNING: found ("+tmpu+") selected ("+(", ".join(uniq))+")"#,vals
		gotwarn = True
	#prepare link with label
	link = "\t"+k+" [label=\""+(",".join(uniq))+"\"];\n"
	fout.write(link)
if(gotwarn): print "To avoid this warning change EMAX in the script."

#select species containing atom base (unique)
uniq = []
for x in species:
	if(x in uniq): continue
	uniq.append(x)
species = uniq


fout.write("{\n")
i = 0
r=10.
for x in species:
	px = r*cos(i*2.*3.1415/(len(species)))
	py = r*sin(i*2.*3.1415/(len(species)))
	fout.write("\""+x+"\" [pos=\""+str(px)+","+str(py)+"\"]\n")
	i += 1
fout.write("}\n")



fout.write("}\n")
fout.close()

#run dot to prepare eps file
fstdout = open(outfile,"w")
try:
	call([prog,tmpfile,"-Teps"],stdout=fstdout)
except:
	print "ERROR: problem with dot!"
	print "You need graphviz to run this script (http://www.graphviz.org/)"
	sys.exit()
os.remove(tmpfile) #delete temporary file


print "Output written in "+outfile
print "DONE, bye!"



