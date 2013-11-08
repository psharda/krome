#	Makefile forKromeChemistry module 
#             D. Seifried 2013

KromeChemistry += KromeChemistry_data.o KromeChemistry_init.o KromeChemistry.o pchem_mapNetworkToSpecies.o \
             krome_all.o krome_user_commons.o opkdmain.o \
             opkda1.o opkda2.o

KromeChemistry_init.o : KromeChemistry_data.o krome_all.o
KromeChemistry.o : KromeChemistry_data.o krome_all.o
krome_all.o : krome_user_commons.o opkdmain.o opkda1.o opkda2.o
opkdmain.o : opkda1.o opkda2.o
opkda1.o : opkda2.o

