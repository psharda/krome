
!###################################################

program test
  use krome_main !use krome
  use krome_user !use utilities
  use krome_constants
  implicit none
  integer,parameter::nsp=5 !number of species
  real*8::Tgas,dt,x(nsp),spy,volume,t,xold(nsp)

  call krome_init() !init krome
  volume = 1d3 !cm3
  x(:) = 1d-3 * N_avogadro / volume !mol/cm3

  Tgas = 4d3 !gas temperature (K)
  dt = 1d-10 !initial time-step (s)
  t = 0.d0 !time (s)
  print *,"solving..."
  do 
     xold(:) = x(:)
     !call the solver
     call krome(x(:), Tgas, dt)
     dt = dt * 1.5d0
     t = t + dt
     write(66,'(99E12.3e3)') t, x(:) / N_avogadro * volume
     if(sum((xold(:)-x(:))**2)<1d-20) exit
  end do

  print *,"done!"
  print *,"Test OK!"
  print *,"Type"
  print *,"gnuplot> load 'plot.gps'"
  print *,"in gnuplot for graphical results."
   
end program test

