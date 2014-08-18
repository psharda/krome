/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS WITH KROME)
/
/  KROME DEVELOPERS, 2014
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
#include "CosmologyParameters.h"
#include "Gadget.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
double ReturnWallTime();
extern "C" void FORTRAN_NAME(krome_driver)(
	float *d, float *e, float *ge, float *u, float *v, float *w, 
        #KROME_args int *in, int *jn, int *kn,
	hydro_method *imethod,
        int *idual, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, 
	float *dt, float *aye,  
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *gamma, float *fh, float *dtoh);


int grid::SolveRateAndCoolEquations(int RTCoupledSolverIntermediateStep)
{
  /* Return if this doesn't concern us. */
  if (!(MultiSpecies && RadiativeCooling)) return SUCCESS;

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
      #KROME_num

  FLOAT a = 1.0, dadt;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. New routine from KROME */

  if (MultiSpecies)
    if (IdentifySpeciesFieldsKrome(
  #KROME_identify
  ) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }


  /* Compute size of the current grid. */

  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* Get easy to handle pointers for each variable. */

  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  /* Compute total gas energy if using MHD */
  if (HydroMethod == MHD_RK) {
    totalenergy = new float[size];
    float B2;
    for (int n=0; n<size; n++) {
      B2 = pow(BaryonField[B1Num][n],2) + pow(BaryonField[B2Num][n],2) + pow(BaryonField[B3Num][n],2);
      totalenergy[n] = BaryonField[TENum][n] - 0.5*B2/BaryonField[DensNum][n];
    }
  }
  else {
    totalenergy = BaryonField[TENum];
  }


  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
            ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  float afloat = float(a);

  float dtCool = dtFixed;

  /* Call the fortran routine to solve cooling equations. */


  FORTRAN_NAME(krome_driver)(
    density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
    #KROME_baryon
    GridDimension, GridDimension+1, GridDimension+2, 
    &HydroMethod, 
    &DualEnergyFormalism,
    &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
    GridEndIndex, GridEndIndex+1, GridEndIndex+2,
    &dtCool, &afloat, 
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
    &Gamma,
    &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio);

  if (HydroMethod == MHD_RK) {
    float B2, v2;
    for (int n = 0; n < size; n++) {
      B2 = pow(BaryonField[B1Num][n],2) + pow(BaryonField[B2Num][n],2) + pow(BaryonField[B3Num][n],2);

      /* Always trust gas energy in cooling routine */
      if (DualEnergyFormalism) {

	v2 = pow(BaryonField[Vel1Num][n],2) + 
	  pow(BaryonField[Vel2Num][n],2) + pow(BaryonField[Vel3Num][n],2);
	BaryonField[TENum][n] = gasenergy[n] + 0.5*v2 + 0.5*B2/BaryonField[DensNum][n];
      }
      else {
	BaryonField[TENum][n] = totalenergy[n] + 0.5*B2/BaryonField[DensNum][n];
      }
      
    }
    
    delete totalenergy;
  }

  return SUCCESS;

}
