#	Makefile for PrimordialChemistry module

PrimordialChemistry += PrimordialChemistry_data.o pchem_mapNetworkToSpecies.o \
             krome_all.o krome_user_commons.o opkdmain.o \
             opkda1.o opkda2.o

PrimordialChemistry_init.o : PrimordialChemistry_data.o krome_all.o
PrimordialChemistry.o : PrimordialChemistry_data.o krome_all.o
krome_all.o : krome_user_commons.o opkdmain.o opkda1.o opkda2.o
opkdmain.o : opkda1.o opkda2.o
opkda1.o : opkda2.o

