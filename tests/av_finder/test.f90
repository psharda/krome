program main
  use krome_main
  use krome_user
  implicit none
  real*8::x(krome_nmols), Av, G0, d2g

  x(:) = 0d0
  x(krome_idx_H2) = 1d0
  d2g = 1d-2

  call krome_init()

  call krome_set_photoBin_draineLin(3d0, 15d0)

  call krome_load_average_kabs()

  call krome_find_G0_Av(G0, Av, x(:), d2g)

  print *, G0, Av

  print *,"done!"

end program main
