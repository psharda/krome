#include <stdio.h>

#include "krome.h"
#include "krome_user.h"

int main() {
  double Tgas, dt, spy, t, ntot, x_60Fe, Fe60_xi;
  double x[krome_nmols];
  int nstep, i;
  FILE * fp;

  spy = krome_seconds_per_year; // seconds per year

  krome_init();

  /* ----==== Isotope decay off ====---- */
  for(i=0;i < krome_nmols;i++) {
    x[i] = 1e-20; // default abundances
  }
  ntot = 1.0e4; // cm^-3
  x[krome_idx_Hj] = ntot; // H+ initial abundance

  // set iron abundance
  krome_scale_z(x, 0.0);
  x[krome_idx_FEj] = x[krome_idx_FE];
  x[krome_idx_FE] = 0.0;
  x[krome_idx_E] = krome_get_electrons(x);

  // set 60Fe abundance
  krome_set_user_tauh(1.5e6 * spy);
  x_60Fe = 0.0;
  printf("Fe-60 abundance = %E\n", x_60Fe);
  x[krome_idx_60FE] = x_60Fe * x[krome_idx_FEj];
  // set rate for ionisation by isotope decay
  Fe60_xi = 0.0;
  printf("Isotope decay off.\n");
  krome_set_user_xi(Fe60_xi);
  // set heating from isotope decay
  krome_set_user_wergs(36.0 * krome_eV_to_erg);

  Tgas = 1.0e3;

  t = 0.0;
  dt = 1.0e-2 * spy;
  nstep = 0;
  fp = fopen("fort.66","w");
  if (fp == NULL) {
    printf("Cannot open outputfile!\n");
    return 1;
  }
  while (t <= 1.0e8 * spy) {
    dt = dt * 1.1;
    t = t + dt;
    if (nstep % 10 == 0) {
      printf("nstep = %d\n",nstep);
    }

    // call KROME
    krome(x, &Tgas, &dt);
    nstep++;

    // output results to file
    fprintf(fp, "%15.8E%15.8E", t/spy, Tgas);
    for (i = 0;i < krome_nmols;i++){
      fprintf(fp, "%15.8E", x[i]/ntot);
    }
    fprintf(fp, "\n");

    krome_popcool_dump(t/spy, 71);
  }

  printf("Finished. Number of steps = %d\n", nstep);

  /* ----==== Isotope decay on ====---- */
  for(i=0;i < krome_nmols;i++) {
    x[i] = 1e-20; // default abundances
  }
  ntot = 1.0e4; // cm^-3
  x[krome_idx_Hj] = ntot; // H+ initial abundance

  // set iron abundance
  krome_scale_z(x, 0.0);
  x[krome_idx_FEj] = x[krome_idx_FE];
  x[krome_idx_FE] = 0.0;
  x[krome_idx_E] = krome_get_electrons(x);

  // set 60Fe abundance
  krome_set_user_tauh(1.5e6 * spy);
  x_60Fe = 1.0e-6;
  printf("Fe-60 abundance = %E\n", x_60Fe);
  x[krome_idx_60FE] = x_60Fe * x[krome_idx_FEj];
  // set rate for ionisation by isotope decay
  Fe60_xi = 1.0e-10;
  printf("Isotope decay on.\n");
  krome_set_user_xi(Fe60_xi);
  // set heating from isotope decay
  krome_set_user_wergs(36.0 * krome_eV_to_erg);

  Tgas = 1.0e3;

  t = 0.0;
  dt = 1.0e-2 * spy;
  nstep = 0;
  fp = fopen("fort.67","w");
  if (fp == NULL) {
    printf("Cannot open outputfile!\n");
    return 1;
  }
  while (t <= 1.0e8 * spy) {
    dt = dt * 1.1;
    t = t + dt;
    if (nstep % 10 == 0) {
      printf("nstep = %d\n",nstep);
    }

    // call KROME
    krome(x, &Tgas, &dt);
    nstep++;

    // output results to file
    fprintf(fp, "%15.8E%15.8E", t/spy, Tgas);
    for (i = 0;i < krome_nmols;i++){
      fprintf(fp, "%15.8E", x[i]/ntot);
    }
    fprintf(fp, "\n");

    krome_popcool_dump(t/spy, 72);
  }
  fclose(fp);

  printf("Finished. Number of steps = %d\n", nstep);

  return 0;
}
