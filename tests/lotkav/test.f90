
!###################################################
! This test is standar Lotka-Volterra model

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  integer,parameter::nsp=krome_nmols !number of species (common)
  real*8::Tgas,dt,x(nsp),t

  call krome_init() !init krome (mandatory)

  x(:) = 1d-20 !default abundances
  x(krome_idx_x) = 5d0 !initial prey
  x(krome_idx_y) = 2d0 !initial predator

  Tgas = 1d3 !gas temperature (dummy)

  dt = 1d-2 !time-step (s)
  t = 0d0 !time (s)
  write(66,'(a)') "#time predators preys"
  do
     call krome(x(:),  Tgas, dt) !call KROME
     write(66,'(99E17.8e3)') t,x(:)
     t = t + dt
     if(t>1d2) exit
  end do
  print *,"done!"
  print *,"to plot type in gnuplot"
  print *,"load 'plot.gps'"
  print *,"have a nice day"

end program test

