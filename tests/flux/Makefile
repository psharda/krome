fc = ifort
#fc = gfortran

lib = #-llapack
ifeq ($(fc),ifort)
	std = -check all -traceback -fpe0  -ftz -ftrapuv -warn all -u
	nowarn = -nowarn
	sswitch = $(std)
	fswitch = $(std) -O3 -ipo -ip -unroll 
	hswitch =  -O3 -ipo -ip -unroll
	pswitch = $(hswitch) -pg
	nswitch = $(std) -O0
	switch =  $(hswitch)
else
	switch += -ffpe-trap=invalid,zero,overflow -g
	switch += -fbacktrace -O3
	nowarn = 
endif



prj = krome

all:
	$(fc) -c opkda2.f $(switch) $(nowarn)
	$(fc) -c opkda1.f $(switch) $(nowarn)
	$(fc) -c opkdmain.f $(switch) $(nowarn)
	$(fc) -c krome_commons.f90 $(switch) 
	$(fc) -c krome_constants.f90 $(switch)
	$(fc) -c krome_user_commons.f90 $(switch)
	$(fc) -c krome_subs.f90 -O0 #$(switch)
	$(fc) -c krome_photo.f90 $(switch)
	$(fc) -c krome_tabs.f90 $(switch)
	$(fc) -c krome_dust.f90 $(switch)	
	$(fc) -c krome_cooling.f90 $(switch)
	$(fc) -c krome_heating.f90 $(switch)	
	$(fc) -c krome_ode.f90 $(switch) 
	$(fc) -c krome_user.f90 $(switch)
	$(fc) -c krome_reduction.f90 $(switch)
	$(fc) -c krome.f90 $(switch) 
	$(fc) -c test.f90 $(switch) 

	$(fc) -o $(prj) *.o $(switch) $(lib)

clean: 
	rm -f *.o $(prj) *~ *.mod *_genmod.f90 *.i90
