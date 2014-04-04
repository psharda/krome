#!/usr/bin/python
#!python

#THIS SCRIPT IS THE DRIVER FOR THE CLASS tabcool
# which is intended to build custom cooling tables
#written by Tommaso Grassi and the KROME team (Apr3,2014)

import tabcool

#filename
fname = "coolCO.dat"

#amount of colliders
xH2 = 1e0
xcoll = {"H2pa": 0.25e0*xH2, "H2or": 0.75*xH2}

#temperature steps
imax = 30


#create the object
tabc = tabcool.krome_tabcool(fname,True)

#init some Tgas
someTgas = [2e1, 2e2, 3e2, 1e3, 2e3]

#loop on Tgas to retrieve cooling in erg/s/cm3
for Tgas in someTgas:
	absdvdz = 1e99 #velocity gradient cm/s/cm
	thin = tabc.get_cool(Tgas,xcoll,absdvdz)

	absdvdz = 1e-11
	thick = tabc.get_cool(Tgas,xcoll,absdvdz)
	
	print Tgas, thin, thick
