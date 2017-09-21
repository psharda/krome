####################
# This example shows how to use the classes of DOCMAKE for
# other purposes than documentation
# it loads a chemical network and perform some action
####################

import sys
#move to the python source folder
sys.path.insert(0, "./src_py/")

import network,options

#load options
myOptions = options.options("options.opt")

#load network
ntw = network.network(myOptions)

#take first reaction
rea = ntw.reactions[0]

#print verbatim
print rea.verbatim

#print reactant names
print "reactants:", [x.name for x in rea.reactants]

#print reactant names
print "reactant mass, g: ", [x.mass for x in rea.reactants]

#produce rate plot
rea.plotRate(myOptions, "plot.png")

#get a species
sp = ntw.getSpeciesByName("H3+")

#print some information
print "name: ", sp.name
print "html: ", sp.nameHtml
print "LaTex: ", sp.nameLatex
print "charge: ", sp.charge
print "enthalpy @ 298.15 K: ", sp.getEnthalpy(ntw.thermochemicalData), "kJ/mol"

#plot cross section to file
sp.plotXsec("xsec.png")

#print photodissociation rate
print "photorate: ", sp.getPhotoRate("leiden", "photodissociation", "BB@1e4K"), "1/s"





