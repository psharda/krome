! test for the Av G0 finder
program main
  use krome_main
  use krome_user
  implicit none
  real*8::x(krome_nmols), Av, G0, d2g

  ! init krome
  call krome_init()

  ! init some chemistry
  x(:) = 0d0
  x(krome_idx_H2) = 1d0
  ! dust to gas mass ratio
  d2g = 1d-2

  ! set some radiation
  call krome_set_photoBin_draineLin(3d0, 15d0)

  ! load absorption coefficients and average over binning
  call krome_load_average_kabs()

  ! find G0 and Av
  call krome_find_G0_Av(G0, Av, x(:), d2g)

  print *, G0, Av

  print *,"done!"

end program main
