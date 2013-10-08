!#################################################################
!This test consists of a gas with an initial population of dust,
! both silicon and carbon based, with 10 bins each.
!The gas-to-dust ratio is 10^-5 and the dust distribution
! follows an MRN power-law.
!#################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_dust
  use krome_constants

  integer,parameter::nx=krome_nmols,nd=10*2
  real*8::x(nx),Tgas,t,dt,spy,xH,tend,vgas,xi(nx)
  real*8::xdust(nd),adust(nd),dust_to_gas(2),xdusti(nd)
  integer::i

  spy = 365.*24.*3600. !seconds per year
  Tgas = 1d1 !gas temperature (K)
  xH = 1d0 !Hydrogen density

  !total abundance of dust (C, Si)
  !ndust(:) = (/1d-3, 1d-4/) !cm-3
  dust_to_gas(:) = 1d-5

  !initialize krome
  call krome_init()

  x(:) = 1.d-20 !default species abundance (cm-3)

  !number densities (cm-3)
  x(KROME_idx_H)     = 1d0 * xH  !H
  x(KROME_idx_Hj)    = 1.d-4 * xH   !H+
  x(KROME_idx_H2)    = 1.d-5 *xH   !H2
  x(KROME_idx_HD)    = 1.d-8 *xH   !HD
  x(KROME_idx_C)     = 1.d-6 *xH    !C
  x(KROME_idx_Si)    = 1.d-6 *xH    !Si

  !compute electrons (globally neutral)
  x(KROME_idx_E) = krome_get_electrons(x(:))

  !initialize dust (xdust=amount per bin, adust=bin size,
  ! ndust=total dust abundance per type, )
  call krome_init_dust(xdust(:), adust(:), dust_to_gas(:),x(:))
  xdusti(:) = xdust(:) !store initial dust amount
  xi(:) = x(:) !store the initial amount of species
  dt = 1d-4*spy !time-step (s)
  t = 0.d0 !initial time (s)
  tend = 1d8*spy !end time (s)
  !write inital values
  write(66,'(999E12.3e3)') t/spy,Tgas,x(:)
  !loop over time-steps
  do
     Tgas = (1d6-1d1) * (t/tend) + 1d1 !Tgas increas with time

     !Tgas = (1d5-1d2) /(1.d0+exp(-25.*(log10(t/spy+1d0)-7.5))) + 1d2
     print '(a10,E11.3,a10,E11.3,a3)',"time:",t/spy,"yr, Tgas:",Tgas,"K"

     call krome(x(:),Tgas,dt,xdust(:)) !###call KROME###

     t = t + dt !increase time
     dt = max(1d2,t/3.d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,Tgas,x(:)/xi(:) !dump species
     !dump dust
     write(77,'(999E12.3e3)') t/spy,adust(nd),Tgas,xdust(nd)/xdusti(nd)
     write(78,'(999E12.3e3)') t/spy,adust(nd/2),Tgas,xdust(nd/2)/xdusti(nd/2)
     write(79,'(999E12.3e3)') t/spy,adust(1),Tgas,xdust(1)/xdusti(1)
     write(80,'(999E12.3e3)') t/spy,adust(nd/2+1),Tgas,&
          xdust(nd/2+1)/xdusti(nd/2+1)
     if(t>tend) exit !exit when overshoot 1d8 years
  end do

  print *,"Output chemistry in fort.66"
  print *,"Output C dust in fort.77"
  print *,"Output Si dust in in fort.78"

  print *,"That's all! have a nice day!"

end program test_krome
