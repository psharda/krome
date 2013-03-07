program test_krome

  use krome_main
  use krome_subs
  use krome_user
  use krome_cooling

  integer,parameter::nspec=12
  integer::i
  real*8::x(nspec),Tgas,rhogas


  rhogas = 1.d0 !dummy value

  !initialize densities
  x(KROME_idx_H)     = 0.9225*rhogas    !H
  x(KROME_idx_E)     = 1.0e-4*rhogas    !E
  x(KROME_idx_Hj)    = 1.0e-4*rhogas    !H+
  x(KROME_idx_D)     = 1.0e-20*rhogas   !D
  x(KROME_idx_Dj)    = 1.0e-20*rhogas   !D+
  x(KROME_idx_HE)    = 0.0972*rhogas    !He
  x(KROME_idx_HEj)   = 1.0e-20*rhogas   !He+
  x(KROME_idx_H2j)   = 1.0e-20*rhogas   !H2+
  x(KROME_idx_H2)    = 1.0e-5*rhogas    !H2
  x(KROME_idx_HD)    = 1.0e-8*rhogas    !HD
  x(KROME_idx_Hk)    = 1.0e-20*rhogas   !H-
  x(KROME_idx_HEjj)  = 1.0e-20*rhogas   !He++

  Tgas = 1d2 !default gas temperature (K)

  !list initial conditions in a table
  call get_infos(x(:), Tgas)

  !create cooling plot in KROME_cooling_plot.dat
  call plot_cool(x)

  print *,"Data saved in KROME_cooling_plot.dat"
  print *,"Default gnuplot: load 'plot.gps'"

  print *,"that's all! have a nice day!"

end program test_krome
