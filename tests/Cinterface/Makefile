#default fortran compiler
fc = ifort
cc = icc
clib = -nofor_main

#default libraries
lib = -llapack

#if dgesv is not present non need for LAPACK
GREP = $(shell grep -i 'dgesv' krome_subs.f90)
ifeq ($(GREP),)
	lib =
endif

#flags
switchOPT = -O3 -ipo -ip -unroll
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused
switchPRO = $(switchOPT) -pg -traceback
switchOMP = $(switchOPT) -openmp

#default switch
switch = $(switchOPT)
switchc = $(switchOPT)

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
objs += krome_dust.o
objs += krome_photo.o
objs += krome_tabs.o
objs += krome_cooling.o
objs += krome_heating.o
objs += krome_ode.o
objs += krome_user.o
objs += krome_reduction.o
objs += krome.o

cobjs = krome_header.o

#default target
#note that Fortran is used as the linker
all: 	$(objs) $(cobjs) testC.o testF.o
	$(fc) $(objs) testF.o -o testF $(switch) $(lib)
	$(fc) $(objs) $(cobjs) testC.o -o test $(switchc) $(clib)

#ifort full debug target
debug: switch = $(switchDBG)
debug: switchc = -traceback
debug: nowarn = $(nowarnIfort)
debug: all

#gfortran target
gfortran: fc = gfortran
gfortran: cc = gcc
gfortran: switch = -ffree-line-length-none
gfortran: switchc =
gfortran: clib = -lm -lgfortran
gfortran: nowarn = $(nowarnGNU)
gfortran: all

#gfortran debug target
gfortran_dbg: fc = gfortran
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: nowarn = $(nowarnGNU)
gfortran_dbg: all

cinterface: $(cobjs)
cinterface: cc = gcc
cinterface: switchc =
cinterface: clib = -lm -lgfortran
cinterface: gfortran

cinterface_opt: $(cobjs)
cinterface_opt: fc = gfortran
cinterface_opt: cc = gcc
cinterface_opt: switch = -O3 -ffree-line-length-none
cinterface_opt: switchc = -O3
cinterface_opt: nowarn = $(nowarnGNU)
cinterface_opt: clib = -lm -lgfortran
cinterface_opt: all

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
	$(cc) $(switchc) $(nowarn) -c $^ -o $@
