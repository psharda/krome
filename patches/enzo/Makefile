#Makefile to compile KROME in ENZO 
fc = ifort
#fc = gfortran

ifeq ($(fc),ifort)
	std = -check all -traceback -fpe0  -ftz -ftrapuv -warn all -u
	sswitch = $(std)
	fswitch = $(std) -O3 -ipo -ip -unroll 
	hswitch =  -O3 -ipo -ip -unroll
	pswitch = $(hswitch) -pg
	nswitch = $(std) -O0
	switch =  $(hswitch)
else
	switch += -ffpe-trap=invalid,zero,overflow -g
	switch += -fbacktrace -O3
endif

all:
	$(fc) -c opkda2.F $(switch) -nowarn
	$(fc) -c opkda1.F $(switch) -nowarn
	$(fc) -c opkdmain.F $(switch) -nowarn
	$(fc) -c krome_user_commons.F90 $(switch)
	$(fc) -c krome_all.F90 $(switch) 

clean: 
	rm -f *.o *~ *.mod *_genmod.f90 *.i90
