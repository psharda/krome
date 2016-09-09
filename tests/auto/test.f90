!###################################################
! This is a test file to understand how automatic reactions
! and bin-based photochemistry work. This is a simple
! network with only
! C+ + e -> C
! C + e -> C+ + e + e
! C + photon -> C+ + e
! the evolution is computed for different j21 values with
! constant Tgas and density. see networks/react_auto file.

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  real*8::Tgas,dt,x(krome_nmols),rho,spy,t,j21s(3),j21
  integer::j

  spy = krome_seconds_per_year !use shorter variable for this constant

  j21s(:) = (/1d-3, 1d-1, 1d0/) !list of j21 values

  call krome_init() !init krome (mandatory)

  !header
  write(66,'(a)') "#j21 time "//krome_get_names_header()
  !loop on j21 values
  do j=1,size(j21s)
     !init radiation field between 1d1 and 1d2 using 1/E profile
     call krome_set_photoBin_J21log(1d1, 1d2)

     j21 = j21s(j) !select J21 value
     call krome_photoBin_scale(j21) !scale radiation according to j21
     print *,"running with j21=",j21
     !init abundances
     x(:) = 1d-40 !default abundances
     x(krome_idx_C) = .5d0
     x(krome_idx_Cj) = .5d0
     x(krome_idx_e) = krome_get_electrons(x(:))

     Tgas = 1.2d4 !gas temperature (K)
     dt = spy !initial time-step (s)
     t = 0d0 !start time (s)

     do
        dt = dt * 1.1d0 !increase time-step
        t = t + dt !advance time
        call krome(x(:), Tgas, dt) !call KROME
        write(66,'(99E17.8)') j21, t/spy, x(:) !dump
        if(t>1d8*spy) exit !exit when 1d6 years is reached
     end do
     write(66,*)
  end do

  print *,"done!"
  print *," for a plot in gnuplot type"
  print *," load 'plot.gps'"

end program test

