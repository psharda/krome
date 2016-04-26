/* C header for krome.f90/the krome_main module.
  Variables for which KROME will return an updated value need to be passed
  by reference using an argument pointer, e.g., "*x".
  Arrays as input should also be passed by reference/argument pointer.
  Passing 2-D arrays as arguments use an array of pointers to an array (i.e., "**x").
  To return an array from a function, it must return an argument pointer,
  e.g., "extern double *functionName()".
*/

#IFNDEF KROME_H_
#DEFINE KROME_H_

#IFKROME_useX
extern void krome(double* x, double rhogas, double *Tgas, double *dt);
extern void krome_equilibrium(double* x, double rhogas, double *Tgas);
#ELSEKROME
extern void krome(double* x, double *Tgas, double *dt);
extern void krome_equilibrium(double* x, double *Tgas);
#ENDIFKROME
extern void krome_init();
extern double *krome_get_coe(double* x, double *Tgas);
extern double *krome_get_coeT(double *Tgas);

#ENDIF
