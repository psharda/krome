!######################################################################
!This test checks the chemical evolution of a stellar interior
! following an already-computed density and temperature evolution
! with FRANEC (e.g. see Tognelli+2011), and evolving the chemistry
! alongside these physical conditions.
!######################################################################
program test
  use krome_main
  use krome_user
  use krome_user_commons
  implicit none
  real*8::x(krome_nmols), Tgas, datar(26), datar2(32), tt, data_pre(32)
  real*8::dtin, rho, T_new, rho_new, time
  integer::ios

  call krome_init()

  !read second line as initial conditions (first line header)
  open(33,file="chem.dat",status="old")
  x(:) = 0.d0 !default abundances
  read(33,*) datar(1:8) !read file line

  !initial abundances
  x(krome_idx_H) = datar(2)
  x(krome_idx_2H) = datar(3)
  x(krome_idx_3He) = datar(4)
  x(krome_idx_4He) = datar(5)
  close(33)

  print *,"running..."
  tt = 0d0 !absolute time (s)
  open(44,file="physcond.dat",status="old")
  read(44,*) data_pre(1:3)
  do
     !reading data
     read(44,*,iostat=ios) time, T_new, rho_new
     if(ios.ne.0) exit

     datar2(1:3) = data_pre(1:3)
     dtin = (10**time - 10**datar2(1)) * krome_seconds_per_year
     data_pre(1) = time
     data_pre(2) = T_new
     data_pre(3) = rho_new
     !write(*,*) datar2(1:3), dtin

     Tgas = 1d1**datar2(2)
     rho = 1d1**datar2(3)

     if(dtin > 0.0) then
        call krome(x(:), rho, Tgas, dtin)
     endif
     write(66,'(F12.8,2F10.5,99(2x,1pe13.6),0p)') datar2(1), &
          log10(Tgas), log10(rho), x(krome_idx_H), x(krome_idx_2H), &
          x(krome_idx_3He), x(krome_idx_4He)

  end do
  close(44)

  print *,"done!"
  print *,"type:"
  print *," >load 'plot.gps'"
  print *,"in gnuplot to show the results"

end program test
