/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: KROME developers, based on the original ENZO routine
/  date:       2013
/
************************************************************************/
 
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int grid::IdentifySpeciesFieldsKrome(
 #KROME_identify
){

 #KROME_zero
 
  /* Find Fields for the KROME species */
#KROME_num

  /* Error if any not found. */
  if (#KROME_binary) {
    ENZO_VFAIL(#KROME_vfail1,
	    #KROME_vfail2)
  }
 
  return SUCCESS;
}
