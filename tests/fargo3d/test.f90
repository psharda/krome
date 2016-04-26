program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  integer,parameter::nsp=krome_nmols !number of species (common)
  real*8::Tgas,dt,x(nsp),rho,spy,ntot

  spy = 3.65d2 * 2.4d1 * 3.6d3 !seconds per year

  call krome_init() !init krome (mandatory)

  !total number density (1/cm3)
  ntot = 1d5

  !init abundances
  x(:) = 1d-20 !default abundances
  x(KROME_idx_H)         = ntot           !H
  x(KROME_idx_H2)        = 1.0e-6*ntot    !H2
  x(KROME_idx_Hj)        = 1.0e-4*ntot    !H+
  x(KROME_idx_HE)        = 0.0775*ntot    !He

  x(krome_idx_e) = krome_get_electrons(x(:))

  Tgas = 3d2 !gas temperature (K)
  dt = spy*1d6 !time-step (s)

  print *,x(krome_idx_Hj)
  call krome(x(:), Tgas, dt) !call KROME
  print *,x(krome_idx_Hj)

  print *,"Test OK!"

end program test

