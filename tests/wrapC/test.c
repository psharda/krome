/*
  The aim of the present code is to show a simple example
  on how to call Krome from a C program. Note that in krome.h
  we provide some useful common variables similar to the one
  found in krome_user.f90
 */

#include <stdio.h>
#include "krome.h" //contains useful krome commons
int main() {
  double Tgas,dt,spy;
  double x[krome_nmols],xi[krome_nmols];
  int i;
  //declares external fortran funcion krome inisde krome_main module
  extern void krome_main_mp_krome_(double *, double *Tgas,double *dt);
  
  spy = 365. * 24. * 3600.; //seconds per year
  
  //default value for species
  for(i=0;i<krome_nmols;i++){
    x[i] = 1e-40;
  }
  x[krome_idx_H] = 1e-6; //H cm-3
  x[krome_idx_H2] = 1e-6; //H2 cm-3
  x[krome_idx_Hj] = 1e0; //H+ cm-3
  x[krome_idx_E] = 1e0; //e- cm-3
  dt = 1e8 * spy; //time-step s
  Tgas = 3e2; //gas temperature K

  //copy initial output values
  for(i=0;i<krome_nmols;i++){
    xi[i] = x[i];
  }

  //call krome!
  krome_main_mp_krome_(x, &Tgas, &dt);

  //print output values
  printf("spec: initial -> final\n");
  printf("----------------------\n");
  for(i=0;i<krome_nmols;i++){
    printf("%s: %e -> %e\n",krome_names[i],xi[i],x[i]);
  }

  printf("\n");
  printf("done!\n");
  return 0;
}
