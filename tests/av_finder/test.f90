! test for the Av G0 finder
program main
  use krome_main
  use krome_user
  implicit none
  integer,parameter::ncell=256
  real*8::r, ntot, dr, rold, rmax, rmin
  real*8::x(krome_nmols), Av, G0, d2g, Av_f
  real*8::tau(krome_nPhotoBins), Tgas
  integer::i

  ! init krome
  call krome_init()

  ntot = 1d3
  Tgas = 1d1
  ! init some chemistry
  x(:) = 0d0
  x(krome_idx_H2) = ntot
  ! dust to gas mass ratio
  d2g = 1d-2

  ! set some radiation
  call krome_set_photoBin_draineLin(5d0, 13.6d0)

  call krome_load_opacity_table("kabs_draine_Rv31.dat", "micron")

  ! load absorption coefficients and average over binning
  call krome_load_average_kabs()

  rmin = 1d10
  rmax = 1d20
  rold = 0d0
  do i=1,ncell
     r = 1d1**((i-1)*(log10(rmax)-log10(rmin))/(ncell-1)+log10(rmin))
     dr = r-rold
     ! find G0 and Av
     tau(:) = krome_get_opacity_size_d2g(x(:), Tgas, dr, d2g)
     call krome_photoBin_scale_array(exp(-tau(:)))
     call krome_find_G0_Av(G0, Av_f, x(:), d2g)

     Av = 2. * ntot * r / 1.8d21 ! there are 2 Hnuclei per H2 molecule(!)
     write(66, '(99E17.8e3)') r, G0, Av_f, Av
     rold = r
  end do

  print *,"done!"

end program main
