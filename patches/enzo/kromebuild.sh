#!/bin/sh

#script to compile KROME in ENZO 
fc=ifort

std="-check all -traceback -fpe0  -ftz -ftrapuv -warn all -u"
hswitch="-O3 -ipo -ip -unroll"
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

echo "removing mod files" 
rm *.mod
echo "everything done!"