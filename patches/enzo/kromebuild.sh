#!/bin/sh

#First copy the new files into ../src/enzo

echo "copying Grid_SolveRateAndCoolEquations.C"

cp Grid_SolveRateAndCoolEquations.C ../

echo "copying Grid_IdentifySpeciesFieldsKrome.C"

cp Grid_IdentifySpeciesFieldsKrome.C ../

echo "copying InitializeRateData.C"

cp InitializeRateData.C ../

#script to compile KROME in ENZO
#WARNING: THIS IS NOT A MAKEFILE!
fc=ifort

std="-check all -traceback -fpe0  -ftz -ftrapuv -warn all -u"
hswitch="-O3"
switch=$hswitch

echo "build using $fc -c $switch"

echo "building opkda2.F"
$fc -c opkda2.F $switch -nowarn
echo "building opkda1.F"
$fc -c opkda1.F $switch -nowarn
echo "building opkdamain.F"
$fc -c opkdmain.F $switch -nowarn
echo "building krome_user_commons.F90"
$fc -c krome_user_commons.F90 $switch
echo "building krome_all.F90"
$fc -c krome_all.F90 $switch

echo "everything done!"
