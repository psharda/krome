#default fortran compiler
fc = ifort
cc = icc

#default libraries
lib = -llapack

#if dgesv is not present non need for LAPACK
GREP = $(shell grep -i 'dgesv' krome_subs.f90)
ifeq ($(GREP),)
	lib =
endif

#flags
switchOPT = -O3 -ipo -ip -unroll -fPIC
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused
switchPRO = $(switchOPT) -pg -traceback
switchOMP = $(switchOPT) -openmp

#default switch
switch = $(switchOPT)

#no warning switch (for DLSODES)
nowarnIfort = -nowarn
nowarnGNU = -w

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += krome_commons.o
objs += krome_constants.o
objs += krome_user_commons.o
objs += krome_fit.o
objs += krome_getphys.o
objs += krome_gadiab.o
objs += krome_grfuncs.o
objs += krome_phfuncs.o
objs += krome_subs.o
objs += krome_stars.o
objs += krome_tabs.o
objs += krome_dust.o
objs += krome_photo.o
objs += krome_cooling.o
objs += krome_heating.o
objs += krome_ode.o
objs += krome_user.o
objs += krome_reduction.o
objs += krome.o

#default target
all: 	$(objs) testF.o
	$(fc) $(objs) testF.o -o testF $(switch) $(lib)
	$(fc) -o libkrome.so $(objs) $(switch) $(lib) -shared

#ifort full debug target
debug: switch = $(switchDBG) -fPIC
debug: nowarn = $(nowarnIfort)
debug: all

#gfortran target
gfortran: fc = gfortran
gfortran: cc = gcc
gfortran: switch = -ffree-line-length-none -fPIC
gfortran: nowarn = $(nowarnGNU)
gfortran: clib = -lm -lgfortran
gfortran: all

#gfortran debug target
gfortran_dbg: fc = gfortran
gfortran_dbg: cc = gcc
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: switchc = switch
gfortran_dbg: nowarn = $(nowarnGNU)
gfortran_dbg: clib = -lm -lgfortran
gfortran_dbg: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@

#rule for c
%.o:%.c
	$(cc) $(switch) $(nowarn) -c $^ -o $@
