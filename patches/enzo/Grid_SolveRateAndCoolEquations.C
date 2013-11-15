/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS)
/
/  RE-ARRANGED BY KROME DEVELOPERS, 2013
/  FROM THE ORIGINAL ENZO ROUTINE
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

/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */

extern int RadiationFieldRecomputeMetalRates;

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);
double ReturnWallTime();
extern "C" void FORTRAN_NAME(krome_driver)(
	float *d, float *e, float *ge, float *u, float *v, float *w, 
 #KROME_args int *in, int *jn, int *kn,
	hydro_method *imethod,
        int *idual, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, 
	float *dt, float *aye, float *temstart, 
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

  this->DebugCheck("SolveRadiativeCooling");

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

  /* Find photo-ionization fields */
  /* NOT NEEDED BY KROME_DRIVER*/
/*  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);
*/

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

  /* Metal cooling codes. */

  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE) {
    if (debug)
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
  }

  /* If both metal fields (Pop I/II and III) exist, create a field
     that contains their sum */

  float *MetalPointer;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  /* Calculate the rates due to the radiation field. */

  if (!GadgetEquilibriumCooling) {
    if (RadiationFieldCalculateRates(Time+0.5*dtFixed) == FAIL) {
        ENZO_FAIL("Error in RadiationFieldCalculateRates.");
    }
  }

  float dtCool = dtFixed;

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  FORTRAN_NAME(krome_driver)(
    density, totalenergy, gasenergy, velocity1, velocity2, velocity3,
 #KROME_baryon
    MetalPointer, GridDimension, GridDimension+1, GridDimension+2, 
    &HydroMethod, 
    &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent,
    &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
    GridEndIndex, GridEndIndex+1, GridEndIndex+2,
    &dtCool, &afloat, &CoolData.TemperatureStart, 
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
    &Gamma,
    &CoolData.HydrogenFractionByMass, &CoolData.DeuteriumToHydrogenRatio,
    &CoolData.SolarMetalFractionByMass,
    &ierr);

  if (ierr) {
      fprintf(stdout, "GridLeftEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(stdout, "GridRightEdge = %"FSYM" %"FSYM" %"FSYM"\n",
	      GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
      fprintf(stdout, "GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	      GridDimension[0], GridDimension[1], GridDimension[2]);
      ENZO_FAIL("Error in FORTRAN rate/cool solver!\n");
  }

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

  delete [] TotalMetals;

  return SUCCESS;

}
