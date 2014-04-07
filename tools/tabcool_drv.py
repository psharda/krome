#!/usr/bin/python
#!python

#THIS SCRIPT IS THE DRIVER FOR THE CLASS tabcool
# which is intended to build custom cooling tables
#written by Tommaso Grassi and the KROME team (Apr3,2014)

import tabcool,sys


#temperature steps
imax = 30


#create the object
tabc = tabcool.krome_tabcool()
tabc.load_BASECOL("CO_He_STS50_V1.kij","He")
tabc.load_LAMDA("co.dat")


#init some Tgas
someTgas = [2e1, 2e2, 3e2, 1e3, 2e3]
someHe = [0e0, .24e0]
#loop on Tgas to retrieve cooling in erg/s/cm3
for xHe in someHe:
	#amount of colliders
	xH2 = 1e0
	xcoll = {"H2pa": 0.25e0*xH2, "H2or": 0.75*xH2, "He":xHe}

	for Tgas in someTgas:
		thin = tabc.get_cool(Tgas,xcoll)

		absdvdz = 1e-11  #velocity gradient cm/s/cm
		thick = tabc.get_cool(Tgas,xcoll,absdvdz)
		
		Tbb = 3e3
		cbb = tabc.get_cool(Tgas,xcoll,absdvdz,Tbb)
	
		print Tgas, xHe, thin, thick,cbb
