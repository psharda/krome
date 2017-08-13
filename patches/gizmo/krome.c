#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#include "krome_all.h"

#ifdef KROME

#define ENDRUNVAL 91234

//Main cooling routine for krome
void do_cooling()
{
  int i;
  double uold,unew;
  double dtime;
  if(All.TimeStep==0) return;
  for(i=FirstActiveParticle;i>=0;i=NextActiveParticle[i])
  {
    if(P[i].Type !=0) continue;
    if(P[i].Mass <=0) continue;
    uold  = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);
    dtime = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
    dtime /= All.cf_hubble_a;
    unew=CallKrome(uold,SphP[i].Density*All.cf_a3inv,dtime,SphP[i].Ne,i,0);
    SphP[i].InternalEnergy = SphP[i].InternalEnergyPred = unew;
    SphP[i].Pressure = get_pressure(i);
  }
}

void convert(double* x,double rho,int mode,double* out)
{
  if(mode==0) //from fractions to densities
	krome_x2n(x,rho,out);
  else	      //from densites to fractions
	krome_n2x(x,rho,out);
}

double krome_get_mmw(int target)
{
   return krome_get_mu_x(SphP[target].krome_species);
}

int krome_get_cooling_time(double n[],double rho,double temp,double* cooling_time)
{
	double dt=86400*365*1.e4; //10^4 yr
        double new_temp=temp;
	double x[KROME_NUM_SPECIES];
	convert(n,rho,0,x);
	krome_thermo(x,&new_temp,dt);
	double rate=-(new_temp-temp)/dt; //-dT/dt (so that rate>0 when cooling)
        *cooling_time=(rate>0 ? temp/rate : MAX_REAL_NUMBER);
        return 1;
}

//
// 'mode' -- tells the routine what to do
//
//     0 == solve chemistry and assign new abundances
//     1 == calculate and return cooling time
//     2 == calculate and return temperature
//     3 == calculate and return gamma
//
double CallKrome(double u_old, double rho, double dt, double ne_guess, int target, int mode)
{
    int k;
    double returnval = 0.0;

    double temp,cooling_time,energy;
    double density;
    double temperature_units = pow(All.UnitVelocity_in_cm_per_s,2.0)*PROTONMASS/BOLTZMANN / All.cf_afac1;
    double density_units = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
    double dtime = dt * All.UnitTime_in_s / All.HubbleParam;

    double x[KROME_NUM_SPECIES];
    double krome_tiny=1.0e-20;

    //Copy variables
    for(k=0;k<KROME_NUM_SPECIES;k++)
      x[k] = SphP[target].krome_species[k];

    density=rho*density_units; //Convert in g/cm^3

    //Set all the relevant quantities for the current step

    //Calculate temperature
    temp = u_old * temperature_units * (SphP[target].krome_gamma-1) * krome_get_mu_x(x);

    //Enforce the temperature floor
    if(temp < All.MinGasTemp) temp = All.MinGasTemp;

    //Set metallicity
#ifdef METALS
    krome_set_z(P[target].Metallicity[0]/All.SolarAbundances[0]);
    krome_set_dust_to_gas(P[target].Metallicity[0]/All.SolarAbundances[0]);
#else
    krome_set_z(0.);
    krome_set_dust_to_gas(0.);
#endif

    switch(mode) {
        case 0:  //solve chemistry & update values

            krome(x,density,&temp,&dtime);

	    //Enforces the minimum temperature allowed, when the gas is colded
	    if(temp < All.MinGasTemp) temp = All.MinGasTemp;

            // Update variables
	    SphP[target].krome_gamma = krome_get_gamma_x(x,temp);

	    for(k=0;k<KROME_NUM_SPECIES;k++)
	      SphP[target].krome_species[k] = DMAX(x[k],krome_tiny);

	    //Update energy
    	    energy = temp / krome_get_mu_x(x) / (SphP[target].krome_gamma-1) / temperature_units;

	    //Update the electron number density
            convert(SphP[target].krome_species,density,0,x);
            SphP[target].Ne = x[krome_idx_E]/(x[krome_idx_H]+x[krome_idx_Hj]);

            returnval = energy;
            break;

	case 1:  //cooling time
            if(krome_get_cooling_time(x,density,temp,&cooling_time) == 0) {
                fprintf(stderr, "Error in calculate_cooling_time.\n");
                endrun(ENDRUNVAL);
            }
            returnval = cooling_time;
            break;
        case 2:  //calculate temperature
            returnval = temp;
            break;
        case 3:  //calculate gamma
            returnval = krome_get_gamma_x(x,temp);
            break;
    } //end switch
    
    return returnval;
}
	


void InitKrome()
{
  krome_set_mpi_rank(ThisTask+1);

  krome_init();

  //Set floor temperature
  krome_set_tfloor(All.MinGasTemp);

  SetIonizationKrome(1);

}

void init_species()
{
  int i,k;
  double density,unit_temperature,temp;
  double krome_tiny = 1.0e-20;

  double x[KROME_NUM_SPECIES];


  unit_temperature = pow(All.UnitVelocity_in_cm_per_s,2.0) * PROTONMASS/BOLTZMANN;

#ifdef METALS
  krome_set_z(All.InitMetallicityinSolar); //Uses the same metallicity for all particles
  krome_set_dust_to_gas(All.InitMetallicityinSolar); //Uses the same metallicity for all particles
#else
  krome_set_z(0.);
  krome_set_dust_to_gas(0.);
#endif

  for(i=0;i<NumPart;i++) 
    {
      if(P[i].Type!=0) continue;

      density=SphP[i].Density * All.cf_a3inv * All.HubbleParam * All.HubbleParam * All.UnitDensity_in_cgs;

      for(k=0;k<KROME_NUM_SPECIES;k++)
      	SphP[i].krome_species[k]=krome_tiny;

      //---------------------------------
      //Here you have to put the species initialisation (one by one)
      //Example (Neutral primordial gas):
      SphP[i].krome_species[krome_idx_H]  = 0.76;
      SphP[i].krome_species[krome_idx_HE] = 0.24;
      //---------------------------------

      //Prevents zeroes in krome species
      for(k=0;k<KROME_NUM_SPECIES;k++)
	SphP[i].krome_species[k]=DMAX(SphP[i].krome_species[k],krome_tiny);

      SphP[i].krome_gamma = GAMMA;

      //Balance electron mass fraction for consistency
      krome_consistent_x(SphP[i].krome_species);

      if(All.InitGasTemp==0)
       temp=SphP[i].InternalEnergy*unit_temperature * krome_get_mu_x(SphP[i].krome_species) * GAMMA_MINUS1;
      else
       temp=All.InitGasTemp;

      if(temp < All.MinGasTemp) temp = All.MinGasTemp;

      SphP[i].krome_gamma = krome_get_gamma_x(SphP[i].krome_species,temp);

      //Update the internal energy to ensure the right initial temperature
      SphP[i].InternalEnergy = temp/unit_temperature/krome_get_mu_x(SphP[i].krome_species)/(SphP[i].krome_gamma-1);
      SphP[i].InternalEnergyPred = SphP[i].InternalEnergy = DMAX(SphP[i].InternalEnergy,All.MinEgySpec);

      //Get density to compute electronic number density
      convert(SphP[i].krome_species,density,0,x);
      SphP[i].Ne = x[krome_idx_E]/(x[krome_idx_H]+x[krome_idx_Hj]);

    }


}

//The init flag must be used to decide whether the UV flux must be re-initialised or not
void SetIonizationKrome(int init)
{
    int k;
    double redshift = 0.0;
    //Set redshift
    if(All.ComovingIntegrationOn)
	redshift=1./All.cf_atime-1;

    krome_set_zredshift(redshift); //Set to the desired redshift
    krome_set_tcmb(2.725*(1+redshift));

    //Add the photo-chemistry initialisation here

}

#endif  //KROME
