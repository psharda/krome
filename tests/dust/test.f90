!#################################################################
!This test consists of a gas with an initial population of dust,
! both silicon and carbon based, with 10 bins each.
!The gas-to-dust ratio is 10^-5 and the dust distribution
! follows an MRN power-law.
!#################################################################
program test_krome

  use krome_main
  use krome_user

  integer,parameter::nd=krome_ndust
  real*8::x(krome_nmols),Tgas,t,dt,spy,xH,tend,vgas,xi(krome_nmols)
  real*8::xdust(nd),adust(nd),xdusti(nd)
  integer::i

  spy = 365.*24.*3600. !seconds per year
  Tgas = 1d1 !gas temperature (K)
  xH = 1d0 !Hydrogen density

  !initialize krome
  call krome_init()

  x(:) = 1.d-20 !default species abundance (cm-3)

  !number densities (cm-3)
  x(KROME_idx_H)     = 1d0 * xH  !H
  x(KROME_idx_Hj)    = 1.d-4 * xH   !H+
  x(KROME_idx_H2)    = 1.d-5 *xH   !H2
  x(KROME_idx_C)     = 1.d-6 *xH    !C
  x(KROME_idx_Si)    = 1.d-6 *xH    !Si

  !compute electrons (globally neutral)
  x(KROME_idx_E) = krome_get_electrons(x(:))

  call krome_set_dust_distribution()

  !scale dust using dust/gas ratio
  call krome_scale_dust_gas_ratio(1d-5, x(:))

  !store intial dust distribution
  xdusti(:) = krome_get_dust_distribution()

  xi(:) = x(:) !store the initial amount of species
  dt = 1d-4*spy !time-step (s)
  t = 0.d0 !initial time (s)
  tend = 1d8*spy !end time (s)
  !write initial values
  write(66,'(999E12.3e3)') t/spy,Tgas,x(:)
  !loop over time-steps
  do
     Tgas = (1d6-1d1) * (t/tend) + 1d1 !Tgas increas with time

     !Tgas = (1d5-1d2) /(1.d0+exp(-25.*(log10(t/spy+1d0)-7.5))) + 1d2
     print '(a10,E11.3,a10,E11.3,a3)',"time:",t/spy,"yr, Tgas:",Tgas,"K"

     call krome(x(:),Tgas,dt) !###call KROME###

     t = t + dt !increase time
     dt = max(1d2,t/3.d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,Tgas,x(:)/xi(:) !dump species
     !dump dust
     xdust(:) = krome_get_dust_distribution()
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
