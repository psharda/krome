fc = ifort
cc = icc
#fc = gfortran

lib = -llapack
std = -check all -traceback -fpe0  -ftz -ftrapuv -warn all -u
nowarn = -nowarn
sswitch = $(std)
fswitch = $(std) -O3 -ipo -ip -unroll
hswitch =  -O3 -ipo -ip -unroll
pswitch = $(hswitch) -pg
nswitch = $(std) -O0
switch =  $(hswitch)
wrap = -lifcore -limf #wrapper flags

prj = test

all:
	$(fc) -c opkda2.f $(switch) $(nowarn)
	$(fc) -c opkda1.f $(switch) $(nowarn)
	$(fc) -c opkdmain.f $(switch) $(nowarn)
	$(fc) -c krome_commons.f90 $(switch)
	$(fc) -c krome_constants.f90 $(switch)
	$(fc) -c krome_user_commons.f90 $(switch)
	$(fc) -c krome_subs.f90 $(switch)
	$(fc) -c krome_photo.f90 $(switch)
	$(fc) -c krome_tabs.f90 $(switch)
	$(fc) -c krome_dust.f90 $(switch)
	$(fc) -c krome_cooling.f90 $(switch)
	$(fc) -c krome_heating.f90 $(switch)
	$(fc) -c krome_ode.f90 $(switch)
	$(fc) -c krome_user.f90 $(switch)
	$(fc) -c krome_reduction.f90 $(switch)
	$(fc) -c krome.f90 $(switch)
	$(cc) -c test.c $(switch)

	$(cc) -o $(prj) *.o $(switch) $(lib) $(wrap)

debug: all

clean:
	rm -f *.o $(prj) *~ *.mod *_genmod.f90 *.i90
