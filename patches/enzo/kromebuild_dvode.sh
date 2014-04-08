#!/bin/sh

#script to compile KROME in ENZO
#WARNING: THIS IS NOT A MAKEFILE!
fc=ifort

std="-check all -traceback -fpe0  -ftz -ftrapuv -warn all -u"
hswitch="-O3"
switch=$hswitch

echo "build using $fc -c $switch"

echo "building dvode_f90_m.F90"
$fc -c dvode_f90_m.F90 $switch -nowarn
echo "building krome_user_commons.F90"
$fc -c krome_user_commons.F90 $switch
echo "building krome_all.F90"
$fc -c krome_all.F90 $switch

echo "everything done!"
