#!/bin/sh

#FIRST SET THE ENZO MAKEFILES TO WORK WITH KROME

echo "copying and building Make.config.objects"

cp ../Make.config.objects ./

echo "" >> Make.config.objects
echo "OBJS_KROME = \\" >> Make.config.objects
echo "       krome_enzo_patch/opkda2.o \\" >> Make.config.objects
echo "       krome_enzo_patch/opkda1.o \\" >> Make.config.objects
echo "       krome_enzo_patch/opkdmain.o \\" >> Make.config.objects
echo "       krome_enzo_patch/krome_initab.o \\" >> Make.config.objects
echo "       krome_enzo_patch/krome_user_commons.o \\" >> Make.config.objects
echo "       krome_enzo_patch/krome_all.o \\" >> Make.config.objects
echo "       krome_enzo_patch/evaluate_tgas.o \\" >> Make.config.objects
echo "       krome_enzo_patch/krome_driver.o" >> Make.config.objects

echo "moving Make.config.objects back"
mv Make.config.objects ../

######
echo "copying and building Make.config.assemble"

cp ../Make.config.assemble ./

sed '
/$(OBJS_ECUDA_LIB)/ {
a\
\               $(OBJS_KROME) \\
n
}' Make.config.assemble > tmp1


cp tmp1 Make.config.assemble
rm tmp1

echo "moving Make.config.assemble back"

mv Make.config.assemble ../

echo "removing temporary file"

echo "done"

echo "you can now compile enzo as usual"
