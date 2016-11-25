import photorates
from math import pi

#create class
p = photorates.photorates()

#access to data (e.g. energy threshold)
Eth_C = p.data["C"]["E_th"]

#load a flux from file (format: energy/eV,flux/[eV/cm2])
p.loadFlux("../../data/black87_incident_eV.dat",Emin=13.6e0,Emax=100.)

#divide by 4pi to have flux per solid angle unit
p.scaleFlux(1e0/4e0/pi)

#atom list
atoms = ["H","He","He+","C","O"]

#loop on atoms
print "Photorates (1/s) with Black 1987 flux in range E=13.6,100 eV"
for atom in atoms:
	print atom + "\t" + str(p.calcRate(atom))


#compute rate, black-body @ 1e3K
p.setBB(3e4)
print "C photorate (1/s) with BB @ 3e4K, ", p.calcRate("C")

#compute rate, black-body @ 1e3K
p.setBB(1e5,verbose=False,Emin=13.6,Emax=1e2)
print "O photorate (1/s) with BB @ 1e5K, ", p.calcRate("O")

#compute rate, Draine flux
p.setDraine()
print "Si photorate (1/s) with Draine flux, ", p.calcRate("Si")
