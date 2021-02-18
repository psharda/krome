#include <stdio.h>

#include "krome.h"
#include "krome_user.h"

/* NOTE: This file mimics the behavior of test.f90 but can be used to test, e.g., building and using
 * the shared library with the C interface. */

/* ################################################### */
/*  HELLO KROME test evolves the network */
/*   FK1 -> FK2 -> FK3 */
/*  where FK* are fake species */

int main(void) {
  double Tgas, dt, t;
  double x[krome_nmols];
  int nstep, i;
  FILE * fp;

  krome_init(); /* init krome (mandatory) */

  for(i=0;i < krome_nmols;i++) {
    x[i] = 0.0; /* default abundances */
  }
  x[krome_idx_FK1] = 1.0;

  Tgas = 1e2; /* gas temperature, not used (K) */
  dt = 1e-5; /* time-step (arbitrary) */
  t = 0e0; /* time (arbitrary) */

  fp = fopen("fort.66","w");
  fprintf(fp, "#time  FK1 FK2 FK3\n");
  while (t <= 5e0) {
    dt = dt * 1.1e0; /* increase timestep */
    krome(x, &Tgas, &dt); /* call KROME */
    fprintf(fp, "%17.8E", t);
    for (i = 0;i < krome_nmols;i++){
      fprintf(fp, "%17.8E", x[i]);
    }
    fprintf(fp, "\n");

    t = t + dt; /* increase time */
  }

  printf("Test OK!\n");
  printf("In gnuplot type\n");
  printf(" load 'plot.gps'\n");
}
