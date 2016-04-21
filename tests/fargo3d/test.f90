program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  integer,parameter::nsp=krome_nmols !number of species (common)
  real*8::Tgas,dt,x(nsp),rho,spy,xtot

  spy = 3.65d2 * 2.4d1 * 3.6d3 !seconds per year

  call krome_init() !init krome (mandatory)

  !total number density (1/cm3)
  xtot = 1d5

  !init abundances
  x(:) = 1d-20 !default abundances
  x(krome_idx_H) = .5d0*xtot
  x(krome_idx_Hj) = .5d0*xtot
  x(krome_idx_e) = krome_get_electrons(x(:))

  Tgas = 3d2 !gas temperature (K)
  dt = spy !time-step (s)

  print *,x(krome_idx_Hj)
  call krome(x(:), Tgas, dt) !call KROME
  print *,x(krome_idx_Hj)

  print *,"Test OK!"

end program test

