program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_dust

  integer,parameter::nx=14,nd=10*2
  real*8::x(nx),Tgas,t,dt,spy,xH,tend
  real*8::xdust(nd),adust(nd),ndust(2)
  integer::i

  spy = 365.*24.*3600. !seconds per year
  Tgas = 1d1 !gas temperature (K)
  xH = 2d4 !Hydrogen density
  
  !total abundance of dust (C, Si)
  ndust(:) = (/1d-3, 1d-4/) !cm-3
  
  !initialize krome
  call krome_init()
  !initialize dust (xdust=amount per bin, adust=bin size,
  ! ndust=total dust abundance)
  call krome_init_dust(xdust(:), adust(:), ndust(:))

  x(:) = 1.d-20 !default species abundance (cm-3)

  !number densities (cm-3)
  x(KROME_idx_H)     = 0.9225d0  !H
  x(KROME_idx_E)     = 1.0d-4    !E
  x(KROME_idx_Hj)    = 1.0d-4    !H+
  x(KROME_idx_D)     = 1.0d-20   !D
  x(KROME_idx_Dj)    = 1.0d-20   !D+
  x(KROME_idx_HE)    = 0.0972d0  !He
  x(KROME_idx_HEj)   = 1.0d-20   !He+
  x(KROME_idx_H2j)   = 1.0d-20   !H2+
  x(KROME_idx_H2)    = 1.0d-5    !H2
  x(KROME_idx_HD)    = 1.0d-8    !HD
  x(KROME_idx_Hk)    = 1.0d-20   !H-
  x(KROME_idx_HEjj)  = 1.0d-20   !He++
  x(KROME_idx_C)     = 1.d-5     !C
  x(KROME_idx_Si)    = 1.d-6     !Si

  !compute electrons (globally neutral)
  x(KROME_idx_e) = x(KROME_idx_Hj) + x(KROME_idx_Dj) + x(KROME_idx_Hej) &
       + x(KROME_idx_H2j) + 2.d0*x(KROME_idx_Hejj) - x(KROME_idx_Hk)

  dt = 1d2*spy !time-step (s)
  t = 0.d0 !initial time (s)
  tend = 1d8*spy !end time (s)
  !write inital values
  write(66,'(999E12.3e3)') t/spy,Tgas,x(:)
  !loop over time-steps
  do
     Tgas = (1d6-1d1) * (t/tend) + 1d1 !Tgas increas with time
     print '(a10,E11.3,a10,E11.3,a3)',"time:",t/spy,"yr, Tgas:",Tgas,"K"

     call krome(x(:),Tgas,dt,xdust(:)) !###call KROME###

     t = t + dt !increase time
     dt = max(1d2,t/3.d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,Tgas,x(:) !dump species
     !dump dust
     do i=1,nd/2 !half C + half Si
        write(77,'(999E12.3e3)') t/spy,adust(i),xdust(i)
        write(78,'(999E12.3e3)') t/spy,adust(nd/2+i),xdust(nd/2+i)
     end do
     write(77,*)
     write(78,*)
     if(t>tend) exit !exit when overshoot 1d8 years
  end do

  print *,"Output chemistry in fort.66"
  print *,"Output C dust in fort.77"
  print *,"Output Si dust in in fort.78"

  print *,"That's all! have a nice day!"

end program test_krome
