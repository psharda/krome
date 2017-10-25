HYDROOBJ += opkda1.o opkda2.o opkdmain.o krome_user_commons.o krome_all.o transfer_ss_amr.o leiden_shielding.o richings_dissociation_rates.o
ifeq ($(strip $(SELFSHIELDING)),yes)
  FFLAGS += -DSELFSHIELDING
endif

#KROME_useLAPACK

NVAR:=$$(( $(NVAR) + #KROME_NMOLS ))
