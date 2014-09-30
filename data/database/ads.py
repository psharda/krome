#BREAK DATABASE
from math import *
mols= ["H","H2","O","O2","O3","OH","CO","CO2","H2O","HO2","H2O2","HCO","H2CO","CH3O","CH3OH"]

def d90(num):
	if(num==0): return "0d0"
	ll = int(log10(num))
	dex = 1e1**ll
	ss = str(round(num/dex,2))+"d"+str(ll)	
	return ss.replace(".0d","d")

data = dict()
data["H"] = [500,650]
data["H2"] = [300,300]
data["O"] = [1700,1700]
data["O2"] = [1250,900]
data["O3"] = [2100,1800]
data["OH"] = [1360,3500]
data["CO"] = [1100,1300]
data["CO2"] = [2300,2300]
data["H2O"] = [4800,4800]
data["HO2"] = [4000,4300]
data["H2O2"] = [6000,5000]
data["HCO"] = [1100,3100]
data["H2CO"] = [1100,3100]
data["CH3O"] = [1100,3100]
data["CH3OH"] = [1100,3100]

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


for rea2 in data2b:
	verb = (" + ".join(rea2[0]))+" -> "+(" + ".join(rea2[1]))
	print "#2-body rate on surface ("+verb+") from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014"
	print "@type: surf2body"
	print "@reacts: "+(",".join([x+"_dust" for x in rea2[0]]))
	print "@prods: "+(",".join([x+"_dust" for x in rea2[1]]))
	print "@limits:"
	Ea = d90(rea2[2])
	r1,r2 = rea2[0]
	Eice1 = d90(data[r1][1])
	Ebare1 = d90(data[r1][0])
	Eice2 = d90(data[r2][1])
	Ebare2 = d90(data[r2][0])
	print "@rate: dust_2body_rate(m(idx_"+r1+"),m(idx_"+r2+"),"+Ea+",krome_dust_asize2(auto_jdust),n(auto_jdust),ice_fraction(auto_jdust),"+Eice1+","+Eice2+","+Ebare1+","+Ebare2+")"
	print


