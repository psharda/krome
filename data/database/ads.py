#BREAK DATABASE
from math import *
import sys
mols= ["H","H2","O","O2","O3","OH","CO","CO2","H2O","HO2","H2O2","HCO","H2CO","CH3O","CH3OH"]

def d90(num):
	if(num==0): return "0d0"
	if(num<1e0): return str(round(num,3))+"d0"
	ll = int(log10(num))
	dex = 1e1**ll
	ss = str(round(num/dex,2))+"d"+str(ll)	
	return ss.replace(".0d","d")

mp = 1.67262178e-24 #g 
me = 9.1093829e-28 #g
mn = 1.674927351e-24 #g 
mnp = mn+mn+me
mH = mp+me
mO = mnp*16e0
mC = mnp*12e0
data = dict()
data["H"] = [500,650,mH]
data["H2"] = [300,300,2e0*mH]
data["O"] = [1700,1700, mO]
data["O2"] = [1250,900,2e0*mO]
data["O3"] = [2100,1800,3e0*mO]
data["OH"] = [1360,3500,mO+mH]
data["CO"] = [1100,1300,mO+mC]
data["CO2"] = [2300,2300,2e0*mO+mC]
data["H2O"] = [4800,4800,2e0*mH+mO]
data["HO2"] = [4000,4300,2e0*mO+mH]
data["H2O2"] = [6000,5000,2e0*(mO+mH)]
data["HCO"] = [1100,3100,mH+mO+mC]
data["H2CO"] = [1100,3100,2e0*mH+mO+mC]
data["CH3O"] = [1100,3100,3e0*mH+mO+mC]
data["CH3OH"] = [1100,3100,4e0*mH+mO+mC]

mlist = sorted([[k,v[2]] for k,v in data.iteritems()], key=lambda x:x[1])

#icount = 0
#for x in mlist:
#	icount += 1
#	print str(icount)+","+x[0]+","+x[0]+"_dust,auto"
#
#for x in mlist:
#	icount += 1
#	print str(icount)+","+x[0]+"_dust,"+x[0]+",auto"
#


data2b=[]
data2b.append([["H","H"],["H2"],0])
data2b.append([["H","O"],["OH"],0])
data2b.append([["H","OH"],["H2O"],0])
data2b.append([["H","O3"],["OH","O2"],480])
data2b.append([["H","H2O2"],["OH","H2O"],1000])
data2b.append([["H","O2"],["HO2"],200])
data2b.append([["H","H2O"],["OH","H2"],9600])
data2b.append([["H","HO2"],["OH","OH"],0])
data2b.append([["H","CO"],["HCO"],600])
data2b.append([["H","HCO"],["H2CO"],0])
data2b.append([["H","H2CO"],["CH3O"],400])
data2b.append([["H","CH3O"],["CH3OH"],0])
data2b.append([["H","HCO"],["CO","H2"],400])
data2b.append([["H","H2CO"],["HCO","H2"],2250])
data2b.append([["H","CH3O"],["H2CO","H2"],150])
data2b.append([["H","CH3OH"],["CH3O","H2"],3000])
data2b.append([["H","CO2"],["CO","OH"],10000])

data2b.append([["O","O"],["O2"],0])
data2b.append([["O","O2"],["O3"],0])
data2b.append([["O","O3"],["O2","O2"],2300])
data2b.append([["O","H2"],["H","OH"],0])
data2b.append([["O","HO2"],["O2","OH"],4640])
data2b.append([["O","OH"],["O2","H"],0])
data2b.append([["O","CO"],["CO2"],160])
data2b.append([["O","HCO"],["CO2","H"],0])
data2b.append([["O","H2CO"],["CO2","H2"],300])

data2b.append([["OH","OH"],["H2O2"],0])
data2b.append([["OH","H2"],["H","H2O"],2100])
data2b.append([["OH","CO"],["CO2","H"],600])
data2b.append([["OH","HCO"],["CO2","H2"],0])
data2b.append([["HO2","H2"],["H","H2O2"],5000])

data2b = sorted(data2b,key=lambda x:len(x[1]))

#icount = 2*len(data)
#fmtold = 0
#for rea2 in data2b:
#	icount += 1
#	rr1 = [x+"_dust" for x in rea2[0]]
#	rr2 = [x+"_dust" for x in rea2[1]]
#	verb = str(icount)+","+(",".join(rr1))+","+(",".join(rr2))+",auto"
#	fmt = len(rea2[1])
#	if(fmt!=fmtold): print "@format:idx,R,R,"+(",".join(["P"]*fmt))+",rate"
#	fmtold = fmt
#	print verb
#sys.exit()



#@type: adsorption
#@reacts: 
#@prods: H+, E
#@limits: 1.360d+01, 5.000d+04
#@rate: sigma_v96(energy_ev, 4.298d-01, 5.475d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
#for mol in mols:
#	print "#Adsorption rate for "+mol+" from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014"
#	print "@type: adsorption"
#	print "@reacts: "+mol
#	print "@prods: "+mol+"_dust"
#	print "@limits:"
#	print "@rate: dust_adsorption_rate(n(idx_"+mol+"), n(nmols+jdust), imsqrt(idx_"+mol+"),ads_stick(jdust),krome_adust2(jdust),sqrTgas)"
#	print


#for mol in mols:
#	print "#Desorption rate for "+mol+" from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014"
#	print "@type: desorption"
#	print "@reacts: "+mol+"_dust"
#	print "@prods: "+mol
#	print "@limits:"
#	Eice = d90(data[mol][1])
#	Ebare = d90(data[mol][0])
#	print "@rate: dust_desorption_rate(fice(auto_jdust),"+Eice+","+Ebare+")"
#	print


print "@var:[ndust] ice_fraction = dust_ice_fraction_array(krome_dust_asize2(:),n(nmols+1:nmols+ndust),n(idx_H2O_dust_1:idx_H2O_dust_1+nmols))"
print "@var:[nspec] m = get_mass()"
print 
print

for rea2 in data2b:
	verb = (" + ".join(rea2[0]))+" -> "+(" + ".join(rea2[1]))
	print "#2-body rate on surface ("+verb+") from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014"
	print "@type: surf2body"
	print "@reacts: "+(",".join([x+"_dust" for x in rea2[0]]))
	print "@prods: "+(",".join([x+"_dust" for x in rea2[1]]))
	print "@limits:"
	Ea = float(rea2[2]) #K
	r1,r2 = rea2[0]
	Eice1 = d90(data[r1][1])
	Ebare1 = d90(data[r1][0])
	Eice2 = d90(data[r2][1])
	Ebare2 = d90(data[r2][0])
	kboltzmann = 1.3806488e-16 #erg/K
	hplanck = 6.6260755e-27 #erg*s
	aa = 1e-8 #cm
	m1 = data[r1][2]
	m2 = data[r2][2]
	mred = m1*m2/(m1+m2)
	P = exp(-aa/3.1415/hplanck*sqrt(2e0*mred*kboltzmann*Ea))
	#print P
	P = d90(P)
	print "@rate: dust_2body_rate("+P+",krome_dust_asize2(auto_jdust-nmols),n(auto_jdust),ice_fraction(auto_jdust-nmols),"+Eice1+","+Eice2+","+Ebare1+","+Ebare2+\
		",krome_dust_T(auto_jdust-nmols))"
	print


