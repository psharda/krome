!######################################################################
!This test checks the chemical evolution of a stellar interior
! following an already-computed density and temperature evolution
! with FRANEC (e.g. see Tognelli+2011), and evolving the chemistry
! alongside these physical conditions.
!######################################################################
program test
  use krome_main
  use krome_user
  implicit none
  real*8::x(krome_nmols), Tgas, datar(26), datar2(32), tt
  real*8::dtin, rho
  integer::ios

  call krome_init()

  !read second line as initial conditions (first line header)
  open(33,file="chemistry.dat",status="old")
  read(33,*) !skip header
  x(:) = 0.d0 !default abundances
  read(33,*) datar(:) !read file line
  x(krome_idx_H) = datar(3)
  x(krome_idx_3He) = datar(5)
  x(krome_idx_4He) = datar(6)
  x(krome_idx_12C) = datar(7)
  x(krome_idx_14N) = datar(8)
  x(krome_idx_16O) = datar(9)
  x(krome_idx_6Li) = datar(11)
  x(krome_idx_7Li) = datar(12)
  x(krome_idx_9Be) = datar(13)
  x(krome_idx_11B) = datar(14)
  close(33)

  print *,"running..."
  rho_old = 1d0
  tt = 0d0 !absolute time (s)
  open(44,file="phy_params.dat",status="old")
  read(44,*) !skip header
  do
     read(44,*,iostat=ios) datar2(:)
     if(ios.ne.0) exit
     Tgas = 1d1**datar2(5)
     rho = 1d1**datar2(7)
     dtin = (1d1**datar2(2))*krome_seconds_per_year - tt
     call krome(x(:), rhogas, Tgas, dtin)
     write(66,'(99E17.8e3)') tt/krome_seconds_per_year, rho, x(:)
     tt = tt + dtin
  end do
  close(44)

  print *,"done!"
  print *,"type:"
  print *," >load 'plot.gps'"
  print *,"in gnuplot to show the results"

end program test
