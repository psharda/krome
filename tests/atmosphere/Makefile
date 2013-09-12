fc = ifort
#fc = gfortran

#lib = -llapack
ifeq ($(fc),ifort)
	std = -check all -traceback -fpe0  -ftz -ftrapuv -warn all -u
	nowarn = -nowarn
	sswitch = $(std)
	fswitch = $(std) -O3 -ipo -ip -unroll 
	hswitch =  -O3 -ipo -ip -unroll
	pswitch = $(hswitch) -pg
	nswitch = $(std) -O0
	switch =  $(nswitch) -warn nounused
else
	switch += -ffpe-trap=invalid,zero,overflow -g
	switch += -fbacktrace -O3
	nowarn = 
endif



prj = test

all:
	$(fc) -c mods.f90 $(switch) 
	$(fc) -c heatlib.f90 $(switch) 
	$(fc) -c test.f90 $(switch) 

	$(fc) -o $(prj) *.o $(switch) $(lib)

clean: 
	rm -f *.o $(prj) *~ *.mod *_genmod.f90 *.i90
