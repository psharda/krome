#BREAK DATABASE
from math import *
import sys
mols= ["H","H2","O","O2","O3","OH","CO","CO2","H2O","HO2","H2O2","HCO","H2CO","CH3O","CH3OH"]

def d90(num,nfrac=3):
	if(num==0): return "0d0"
	#if(num<1e0): return str(round(num,3))+"d0"
	ll = int(log10(num))
	dex = 1e1**ll
	if(num<1e0):
		zrs = "0"*(nfrac+1-len(str(round(num/dex,nfrac+1)).split(".")[1]))
		ss = str(round(num/dex,nfrac+1)*1e1)+zrs+"d"+str(ll-1)
	else:
		zrs = "0"*(nfrac-len(str(round(num/dex,nfrac)).split(".")[1]))
		ss = str(round(num/dex,nfrac))+zrs+"d"+str(ll)

	return ss #.replace(".0d","d")

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

data2b=[]
data2b.append([["H","H"],["H2"],0,0e0])
data2b.append([["H","O"],["OH"],0,.5])
data2b.append([["H","OH"],["H2O"],0,.9])
data2b.append([["H","O3"],["OH","O2"],480,0e0])
data2b.append([["H","H2O2"],["OH","H2O"],1000,0e0])
data2b.append([["H","O2"],["HO2"],200,0e0])
data2b.append([["H","H2O"],["OH","H2"],9600,0e0])
data2b.append([["H","HO2"],["OH","OH"],0,0e0])
data2b.append([["H","CO"],["HCO"],600,0e0])
data2b.append([["H","HCO"],["H2CO"],0,0e0])
data2b.append([["H","H2CO"],["CH3O"],400,.5])
data2b.append([["H","CH3O"],["CH3OH"],0,.05])
data2b.append([["H","HCO"],["CO","H2"],400,.05])
data2b.append([["H","H2CO"],["HCO","H2"],2250,0e0])
data2b.append([["H","CH3O"],["H2CO","H2"],150,0e0])
data2b.append([["H","CH3OH"],["CH3O","H2"],3000,0e0])
data2b.append([["H","CO2"],["CO","OH"],10000,0e0])

data2b.append([["O","O"],["O2"],0,.6])
data2b.append([["O","O2"],["O3"],0,0e0])
data2b.append([["O","O3"],["O2","O2"],2300,0e0])
data2b.append([["O","H2"],["H","OH"],0,0e0])
data2b.append([["O","HO2"],["O2","OH"],4640,0e0])
data2b.append([["O","OH"],["O2","H"],0,0e0])
data2b.append([["O","CO"],["CO2"],160,0e0])
data2b.append([["O","HCO"],["CO2","H"],0,0e0])
data2b.append([["O","H2CO"],["CO2","H2"],300,0e0])

data2b.append([["OH","OH"],["H2O2"],0,0e0])
data2b.append([["OH","H2"],["H","H2O"],2100,0e0])
data2b.append([["OH","CO"],["CO2","H"],600,0e0])
data2b.append([["OH","HCO"],["CO2","H2"],0,0e0])
data2b.append([["HO2","H2"],["H","H2O2"],5000,0e0])

data2b = sorted(data2b,key=lambda x:len(x[1]))


#****ADSORPTION****
#@type: adsorption
#@reacts: 
#@prods: H+, E
#@limits: 1.360d+01, 5.000d+04
#@rate: sigma_v96(energy_ev, 4.298d-01, 5.475d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
#@rate: dust_adsorption_rate(n(idx_H), n(auto_jdust), imsqrt(idx_H),ads_stick(auto_jdust-nmols),krome_dust_asize2(auto_jdust-nmols),sqrTgas)
print "writing adsorption..."
fout = open("surface_adsorption.dat","w")
fout.write("@var:[nspec] imsqrt = get_imass_sqrt()\n")
fout.write("@var:[ndust] ads_stick = dust_stick_array(Tgas,krome_dust_T)\n")
fout.write("@var: sqrTgas = sqrt(Tgas)\n\n")
for mol in mols:
	fout.write("#Adsorption rate for "+mol+" from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014\n")
	fout.write("@type: adsorption\n")
	fout.write("@reacts: "+mol+"\n")
	fout.write("@prods: "+mol+"_dust\n")
	fout.write("@limits:\n")
	fout.write("@rate: dust_adsorption_rate(xdust(auto_jdust-nmols), imsqrt(idx_"+mol+"),"\
		+"ads_stick(auto_jdust-nmols),krome_dust_asize2(auto_jdust-nmols),sqrTgas)\n")
	fout.write("\n")
fout.close()

#****DESORPTION****
print "writing desorption..."
fout = open("surface_desorption.dat","w")
fout.write("@var:[ndust] dust_inv_phi = dust_get_inv_phi(krome_dust_asize2(:),xdust(:))\n\n")
#fout.write("@var:[ndust] invTdust = 1d0/(krome_dust_T(:)+1d-40)\n")
#fout.write("@var:[2*nspec] Ebareice_exp = get_Ebareice_exp_array(invTdust(:))\n")
fout.write("@var:[2*nspec] Ebareice_exp = dust_Ebareice_exp(:)\n")
fout.write("@var:[nspec] Eice_exp = Ebareice_exp(1:nspec)\n")
fout.write("@var:[nspec] Ebare_exp = Ebareice_exp(nspec+1:2*nspec)\n")
fout.write("@var:[ndust] ice_fraction = dust_ice_fraction_array(dust_inv_phi(:),n(idx_H2O_dust_1:idx_H2O_dust_1+ndust))\n")
for mol in mols:
	fout.write("#Desorption rate for "+mol+" from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014\n")
	fout.write("@type: desorption\n")
	fout.write("@reacts: "+mol+"_dust\n")
	fout.write("@prods: "+mol+"\n")
	fout.write("@limits:\n")
	#Eice = d90(data[mol][1])
	#Ebare = d90(data[mol][0]) #dust_desorption_rate(ice_fraction(auto_jdust-nmols),6.5d2,5d2,krome_dust_T(auto_jdust-nmols))
	Eice = "Eice_exp(idx_"+mol+"_DUST_auto_idx)"
	Ebare = "Ebare_exp(idx_"+mol+"_DUST_auto_idx)"
	fout.write("@rate: dust_desorption_rate(ice_fraction(auto_jdust-nmols),"+Eice+","+Ebare+")\n")
	fout.write("\n")
fout.close()

#****2BODY****
print "writing 2body..."
fout = open("surface_2body.dat","w")
fout.write("@var:[ndust] dust_inv_phi = dust_get_inv_phi(krome_dust_asize2(:),xdust(:))\n\n")
#fout.write("@var:[ndust] invTdust = 1d0/(krome_dust_T(:)+1d-40)\n")
fout.write("@var:[nspec] m = get_mass()\n")
#fout.write("@var:[2*nspec] Ebareice23_exp = get_Ebareice23_exp_array(invTdust(:))\n")
fout.write("@var:[2*nspec] Ebareice23_exp = dust_Ebareice23_exp(:)\n")
fout.write("@var:[nspec] Eice23_exp = Ebareice23_exp(1:nspec)\n")
fout.write("@var:[nspec] Ebare23_exp = Ebareice23_exp(nspec+1:2*nspec)\n")
fout.write("@var:[ndust] ice_fraction = dust_ice_fraction_array(dust_inv_phi(:),n(idx_H2O_dust_1:idx_H2O_dust_1+nmols))\n")
fout.write("\n\n")

for rea2 in data2b:
	verb = (" + ".join(rea2[0]))+" -> "+(" + ".join(rea2[1]))
	fout.write("#2-body rate surface-surface ("+verb+") from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014\n")
	fout.write("@type: surf2body\n")
	fout.write("@reacts: "+(",".join([x+"_dust" for x in rea2[0]]))+"\n")
	fout.write("@prods: "+(",".join([x+"_dust" for x in rea2[1]]))+"\n")
	fout.write("@limits:\n")
	Ea = float(rea2[2]) #K
	delta_bare = float(rea2[3]) #desorpion probability bare
	delta_ice = delta_bare / 5e0 #desorpion probability ice
	r1,r2 = rea2[0]
	Eice1 = d90(data[r1][1])
	Ebare1 = d90(data[r1][0])
	Eice2 = d90(data[r2][1])
	Ebare2 = d90(data[r2][0])
	kboltzmann = 1.3806488e-16 #erg/K
	hplanck = 6.6260755e-27 #erg*s
	hplanck_bar = hplanck/2e0/3.1415 #erg*s
	aa = 1e-8 #cm
	m1 = data[r1][2]
	m2 = data[r2][2]
	mred = m1*m2/(m1+m2)
	#tunnelling probability
	P = exp(-2e0*aa/hplanck_bar*sqrt(2e0*mred*kboltzmann*Ea))
	P = d90(P)
	Eice1 = "Eice_exp(idx_"+r1+"_DUST_auto_idx)"
	Ebare1 = "Ebare_exp(idx_"+r1+"_DUST_auto_idx)"
	Eice2 = "Eice_exp(idx_"+r2+"_DUST_auto_idx)"
	Ebare2 = "Ebare_exp(idx_"+r2+"_DUST_auto_idx)"
	arg_delta_ice = "1d0-"+d90(delta_ice) 
	arg_delta_bare = "1d0-"+d90(delta_bare)
	if(delta_ice==0e0): arg_delta_ice = arg_delta_bare = "1d0" 
	fout.write("@rate: dust_2body_rate("+P+",dust_inv_phi(auto_jdust-nmols),ice_fraction(auto_jdust-nmols),"\
		+Eice1+","+Eice2+","+Ebare1+","+Ebare2\
		+","+arg_delta_ice+","+arg_delta_bare+")\n\n")


	if(delta_ice==0e0): continue
	fout.write("#2-body rate surface-gas ("+verb+") from Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014\n")
	fout.write("@type: surf2body\n")
	fout.write("@reacts: "+(",".join([x+"_dust" for x in rea2[0]]))+"\n")
	fout.write("@prods: "+(",".join([x for x in rea2[1]]))+"\n")
	fout.write("@limits:\n")
	arg_delta_ice = d90(delta_ice) 
	arg_delta_bare = d90(delta_bare)
	fout.write("@rate: dust_2body_rate("+P+",dust_inv_phi(auto_jdust-nmols),ice_fraction(auto_jdust-nmols),"\
		+Eice1+","+Eice2+","+Ebare1+","+Ebare2\
		+","+arg_delta_ice+","+arg_delta_bare+")\n\n")


fout.close()
#*************CHEMISORPTION FUNCTIONS*******************
#eqn. 1+3 from Cazaux+Tielens 2004 (see erratum)
def Tij1(xvar,Bi,Bj,Bij,Z,Tsys):

	hbar = 1.05457266e-27 #erg*s
	kb = 1.380658e-16 #erg/K
	m = 1.6733e-24 #g

	if(xvar>=Bi): return 0e0

	Tij1 = 4e0 *sqrt((xvar-Bij)/xvar) \
		 / ((1e0+sqrt((xvar-Bij)/xvar))**2 \
		 + Bi*Bj*(sinh(Z*sqrt(2e0*m*(Bi-xvar)*kb/hbar**2)))**2 \
		 / (Bi-xvar)/xvar)
	return exp(-xvar/Tsys)*Tij1 / Tsys

#eqn. 2+3 from Cazaux+Tielens 2004 (see erratum)
def Tij2(xvar,Bi,Bj,Bij,Z,Tsys):

	hbar = 1.05457266e-27 #erg*s
	kb = 1.380658e-16 #erg/K
	m = 1.6733e-24 #g

	if(xvar<=Bi): return 0e0

	Tij2 = 4e0 *sqrt((xvar-Bij)/xvar) \
		 / ((1e0+sqrt((xvar-Bij)/xvar))**2 \
		 - Bi*Bj*(sin(Z*sqrt(2e0*m*(xvar-Bi)*kb/hbar**2)))**2 \
		 / (Bi-xvar)/xvar)
	return exp(-xvar/Tsys)*Tij2 / Tsys

#*************CHEMISORPTION*******************
from scipy.integrate import quad
#data from Iqbal+2012 ApJ
Ep = 780 #K
Ec = 1.4e4 #K
Es = 1e2 #K
Esp = 1e2 #K
Esc = 7e3 #K
alow = 2.5e-8 #cm
aup = 2e-8 #cm
nu0 = 1e12 # 1/s

datachemis = []
datachemis.append(["PP",Ep-Esp,Ep-Esp,0e0,aup])
datachemis.append(["CC",Ec-Esc,Ec-Esc,0e0,aup])
datachemis.append(["PC",Ep-Es,Ec-Es,Ep-Ec,alow])
datachemis.append(["CP",Ec-Es,Ep-Es,Ep-Ec,alow])
print "computing chemisorption..."
print "writing rates in ../. (data)..."
fout = open("../surface_chemisorption_rates.dat","w")
fout.write("#Rates for H chemisorption rates. Lines are:\n")
fout.write("#1. process type (P=physisorbed, C=chemisorbed)\n")
fout.write("#2. number of temperature interval (linear)\n")
fout.write("#3. mininum temperature\n")
fout.write("#4. temperature interval\n")
fout.write("#data are rate (1/s)\n\n")

rateChemis = dict()
for datac in datachemis:
	Bi = datac[1]
	Bj = datac[2]
	Bij = datac[3]
	Z = datac[4]
	imax = 300
	Tmin = 1e0
	Tmax = 1e4
	ydata = []
	fout.write(datac[0]+"\n")
	fout.write(str(imax)+"\n")
	fout.write(d90(Tmin)+"\n")
	fout.write(d90((Tmax-Tmin)/imax)+"\n")
	for i in range(imax):
		Tsys = i*(Tmax-Tmin)/(imax-1)+Tmin
		Ptunnel = quad(Tij1, 1e-40, Bi, args=(Bi,Bj,Bij,Z,Tsys),limit=5000,epsabs=1e-40)[0]
		Pdiff = quad(Tij2, Bi, -log(1e-40)*Tsys, args=(Bi,Bj,Bij,Z,Tsys),limit=5000,epsabs=1e-40)[0]
		rate_val = nu0*(Pdiff+Ptunnel)
		ydata.append(rate_val)
		fout.write(d90(rate_val,7)+"\n")
	fout.write("\n")
	rateChemis[datac[0]] = {"rate":ydata, "Tmin":Tmin, "dT":(Tmax-Tmin)/imax,"ndata":imax}

fout.close()

reactChemis = []
reactChemis.append([["H_dust"],["H_c_dust"],["PC"]])
reactChemis.append([["H_c_dust"],["H_dust"],["CP"]])
reactChemis.append([["H_c_dust","H_dust"],["H2_dust"],["CP"]])
reactChemis.append([["H_c_dust","H_dust"],["H2_dust"],["PC"]])
reactChemis.append([["H_c_dust","H_c_dust"],["H2_dust"],["CC"]])
#reactChemis.append([["H_c_dust"],["H"],["CG"]])
print "writing chemisorption..."
fout = open("surface_chemisorption.dat","w")
fout.write("@var:[ndust] rateChem_PC = dust_get_rateChem_PC(krome_dust_T(:))\n")
fout.write("@var:[ndust] rateChem_CP = dust_get_rateChem_CP(krome_dust_T(:))\n")
fout.write("@var:[ndust] rateChem_CC = dust_get_rateChem_CC(krome_dust_T(:))\n")
#fout.write("@var:[ndust] rateChem_CG = dust_get_rateChem_CG(krome_dust_T(:))\n")
fout.write("@var:[ndust] dust_inv_phi = dust_get_inv_phi(krome_dust_asize2(:),xdust(:))\n\n")
for rChem in reactChemis:
	verb = (" + ".join([x.replace("_dust","") for x in rChem[0]]))+" -> "+(" + ".join([x.replace("_dust","") for x in rChem[1]]))
	fout.write("#rate chemisorption ("+verb+") from Cazaux+2004 (err:Cazaux+2010), Iqbal+2010, Hocuk+2014\n")
	fout.write("@type: surfChemisorption\n")
	fout.write("@reacts: "+(",".join(rChem[0]))+"\n")
	fout.write("@prods: "+(",".join(rChem[1]))+"\n")
	fout.write("@limits:\n")
	dust_inv_phi = ""
	if(len(rChem[0])==2): dust_inv_phi = "*dust_inv_phi(auto_jdust-nmols)"
	fout.write("@rate: rateChem_"+rChem[2][0]+"(auto_jdust-nmols)"+dust_inv_phi+"\n\n")

fout.close()


print "done!"
