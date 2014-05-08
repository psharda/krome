#!/usr/bin/python
#!python
from kromelib import *
ys = ["Y","y"]
ns = ["N","n"]
yns = ys + ns
coolList = ["NONE","ATOMIC","H2","HD","DH","DUST","H2GP98","COMPTON","EXPANSION","CIE","CONT","CHEM","DISS","Z"]
coolHelp = ["No cooling","Atomic from Cen 1992","H2 from Glover+2007","HD from Lipovka+2007","Endothermic with thermochemical data","Dust cooling",\
		"H2 from Galli+Palla 1998","Compton","Isothermal expanding gas","Collisional induced","Continuum emission",\
		"Endothermic reactions","Dissociative","Metal"]

heatList = ["COMPRESS", "PHOTO", "CHEM", "DH", "CR", "PHOTOAV", "PHOTODUST"]
heatHelp = ["isothermal copmpressional", "photoheating","exothermic reactions", "exothermic from thermochemical data","cosmic rays",\
	"photoheating due to UV photopumping","photoelectric heating from dust"]

print "********************"
print " option wizard!"
print "********************"
rawc = 0
fullcmd = "./krome"

fname = ""
#READ FILENAME networks/react_COthin
rawc += 1
while(fname==""):
	rwi = raw_input(str(rawc)+". Path of your chemical network [networks/react_COthin]: ").strip()
	if(rwi==""): rwi = "networks/react_COthin"
	if(file_exists(rwi)):
		fname = rwi
	else:
		print rwi +" file not found, try again..."
fullcmd += " -n "+fname

#NUMBER DENSITY
rwi = ""
rawc += 1
while(not(rwi in yns)):
	rwi = raw_input(str(rawc)+'. Use number density (otherwise mass fractions)? [n]: ').strip()
	if(rwi==""): rwi = "n"
	useN = (rwi in ys)
if(useN): fullcmd += " -useN"


#COOLING FILE
rawc += 1
rwi = ""
fcool = ""
while(fcool==""):
	rwi = raw_input(str(rawc)+'. Do you want to extend cooling with an additional cooling file? [data/coolZ.dat]: ').strip()
	if(rwi==""): rwi = "data/coolZ.dat"
	if(file_exists(rwi)):
		fcool = rwi
	else:
		print rwi +" file not found, try again..."

if(fcool!="NONE"): fullcmd += " -coolFile="+fcool

#COOLINGS
rawc += 1
cools = ""
while(cools==""):
	for i in range(len(coolList)):
		print str(i)+")",coolList[i],coolHelp[i]
	rwi = raw_input(str(rawc)+'. Cooling functions (use numbers above comma separated) [0]: ').strip()
	if(rwi==""): rwi = "0"
	cools = ",".join([coolList[int(x)] for x in sorted(rwi.split(","))])
if(cools!="NONE"): fullcmd += " -cooling="+cools


#HEATINGS
rawc += 1
heats = ""
while(heats==""):
	for i in range(len(heatList)):
		print str(i)+")",heatList[i],heatHelp[i]
	rwi = raw_input(str(rawc)+'. Heating functions (use numbers above comma separated) [0]: ').strip()
	if(rwi==""): rwi = "0"
	heats = ",".join([heatList[int(x)] for x in sorted(rwi.split(","))])
if(heats!="NONE"): fullcmd += " -heating="+heats

print
print fullcmd


