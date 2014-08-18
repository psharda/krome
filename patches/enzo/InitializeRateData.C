/***********************************************************************
/
/  INITIALIZE THE MULTI-SPECIES RATES
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  Dan Reynolds, July 2010; added case-B recombination rates
/  modified2:  Britton Smith, October 2010; moved reading/writing of 
/              parameters to Read/WriteParameterFile.
/
/  PURPOSE:
/    For multi-species runs (with cooling), initialize both the
/      CoolData and RateData rate tables.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/* THIS IS WHAT YOU NEED YOU WANT TO USE RATE TABLES FROM KROME */
extern "C" void FORTRAN_NAME(krome_initab)();

int InitializeRateData(FLOAT Time)
{
 
  /* Declarations. */
 
  FLOAT a = 1, dadt;

  /* If using cosmology, compute the expansion factor and get units. */
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt)
	== FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
  }
  float afloat = float(a);
  int ioutput = (( MyProcessorNumber == ROOT_PROCESSOR ) ? 1 : 0);
 
  /* Call FORTRAN routine to prepare tables: KROME 2013 */

    FORTRAN_NAME(krome_initab)();

  return SUCCESS;
}
