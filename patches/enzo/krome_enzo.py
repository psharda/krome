import os
from shutil import copyfile, move

fname_make_config_objects = "Make.config.objects"
fname_make_config_assemble = "Make.config.assemble"
fname_grid_h = "Grid.h"
fname_gridkrome_h = "GridKrome.h"

#fname_make_config_objects
print "patching ",fname_make_config_objects
copyfile("../"+fname_make_config_objects, "./"+fname_make_config_objects)

fh = open(fname_make_config_objects,"a")
fh.write("\n")
fh.write("OBJS_KROME = \\\n")
fh.write("       Grid_IdentifySpeciesFieldsKrome.o\\\n")
fh.write("       krome_enzo_patch/opkda2.o\\\n")
fh.write("       krome_enzo_patch/opkda1.o\\\n")
fh.write("       krome_enzo_patch/opkdmain.o\\\n")
fh.write("       krome_enzo_patch/krome_initab.o\\\n")
fh.write("       krome_enzo_patch/krome_user_commons.o\\\n")
fh.write("       krome_enzo_patch/krome_all.o\\\n")
fh.write("       krome_enzo_patch/evaluate_tgas.o\\\n")
fh.write("       krome_enzo_patch/krome_driver.o\\\n")
fh.close()

move("./"+fname_make_config_objects, "../"+fname_make_config_objects)

#fname_make_config_assemble
print "patching ",fname_make_config_assemble
copyfile("../"+fname_make_config_assemble, "./"+fname_make_config_assemble)

fh = open(fname_make_config_assemble,"rb")
fout = open("tmp1.krome","w")
for row in fh:
	fout.write(row)
	if("$(OBJS_ECUDA_LIB) \\" in row): 
		fout.write("               $(OBJS_KROME) \\\n")
fh.close()
fout.close()
move("tmp1.krome","../"+fname_make_config_assemble)
os.remove("./"+fname_make_config_assemble)

#fname_grid_h
print "patching ",fname_grid_h
fh = open("../"+fname_grid_h,"rb")
fout = open("tmp1.krome","w")
inBlock = False
for row in fh:
	if("int IdentifySpeciesFields(" in row):
		inBlock = True
		fout.write(row)
		continue
	fout.write(row)
	if(inBlock and (";" in row)):
		inBlock = False
		fh2 = open(fname_gridkrome_h,"rb")
		for row2 in fh2:
			fout.write(row2)
		fh2.close()
fh.close()
fout.close()

move("tmp1.krome","../"+fname_grid_h)
print "done!"

	


