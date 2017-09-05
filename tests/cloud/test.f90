!##########################################################
!The following test provides a benchmark for
! molecular clouds environments and it is useful
! to test the capability of the solver to handle
! with very large network.
!It is based on the EA2 model of Wakelam&Herbst 2008.
!##########################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons

  integer,parameter::nx=krome_nmols
  real*8::x(nx),Tgas,t,dt,spy,xH

  spy = 3600. * 24. * 365. !seconds per year
  Tgas = 1d1 !gas temperature (K)
  xH = 2d4 !Hydrogen density

  !user commons for opacity and CR rate
  call krome_set_user_av(1d1) !opacity Av (#)
  call krome_set_user_crate(1.3d-17) !CR rate (1/s)
  call krome_set_user_gas_dust_ratio(7.57d11) !gas/dust

  call krome_init()

  x(:) = 1.d-20
  !initial densities (model EA2 Wakelam+Herbst 2008)
  x(KROME_idx_H2)  = 0.5d0   * xH
  x(KROME_idx_He)  = 9d-2   * xH
  x(KROME_idx_N)   = 7.6d-5  * xH
  x(KROME_idx_O)   = 2.56d-4 * xH
  x(KROME_idx_Cj)  = 1.2d-4  * xH
  x(KROME_idx_Sj)  = 1.5d-5  * xH
  x(KROME_idx_Sij) = 1.7d-6  * xH
  x(KROME_idx_Fej) = 2d-7   * xH
  x(KROME_idx_Naj) = 2d-7   * xH
  x(KROME_idx_Mgj) = 2.4d-6  * xH
  x(KROME_idx_Clj) = 1.8d-7  * xH
  x(KROME_idx_Pj)  = 1.17d-7 * xH
  x(KROME_idx_Fj)  = 1.8d-8  * xH

  !calculate elctrons (neutral cloud)
  x(KROME_idx_e) = krome_get_electrons(x(:))

  !NOTE: here myCoe array is employed to store the
  ! coefficient values, since the temperature is
  ! constant during the model evolution.
  ! myCoe(:) is defined in krome_user_commons
  myCoe(:) = krome_get_coef(Tgas,x(:))

  dt = 1d2*spy !time-step (s)
  t = 0d0 !initial time (s)

  !output header
  write(66,'(a)') "#time "//trim(krome_get_names_header())

  do
     print '(a10,E11.3,a3)',"time:",t/spy,"yr"
     call krome(x(:),Tgas,dt) !call KROME
     t = t + dt !increase time
     dt = max(1d2,t/3d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,x(:)/2d4
     if(t>1d8*spy) exit !exit when overshoot 1d8 years
  end do

  print *,"Output dump in fort.66"
  print *,"Default gnuplot: load 'plot.gps'"

  print *,"That's all! have a nice day!"

end program test_krome
