!################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a metal-enriched cloud.
!The dynamics is described by the Larson-Penston-type
! similar solution and includes cooling and heating processes.
!For additional details look also to Omukai 2005.
!###############################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  integer,parameter::rstep = 500000
  integer,parameter::nzs=4
  integer::i
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt,Tdust
  real*8::ntot,rho,zs(nzs),amin,amax,d2g

  !INITIALIZE KROME PARAMETERS AND DUST
  call krome_init()

  !$omp parallel do schedule(dynamic,1) default(none) &
  !$omp  private(jz,ntot,Tgas,x,dd,i,dd1,rho,tff,dt,dtH,deldd) &
  !$omp  shared(zs,imax,result)

  !INITIAL CONDITIONS
  ntot           = 0.1d0  !total density in 1/cm3
  Tgas           = 3d2    !temperature in kelvin

  !species initialization in 1/cm3
  x(:) = 1d-40

  x(KROME_idx_H)  = ntot          !H
  x(KROME_idx_H2) = 1d-6*ntot    !H2
  x(KROME_idx_Hj) = 1d-4*ntot    !H+
  x(KROME_idx_He) = 0.0775d0*ntot !He

  !set cosmic rays
  call krome_set_user_crate(1.3d-17)
  !call krome_set_user_Tdust(1d1)

  !rescale metallicity for neutral metals (C,Fe,Si,O)
  call krome_scale_Z(x(:), 0d0)

  x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
  x(krome_idx_C)  = 1d-40

  x(krome_idx_e) = krome_get_electrons(x(:))

  dd = ntot
  print '(2a5,2a11)',"Zstep","step","n(cm-3)","Tgas(K)"

  !loop over the hydro time-step
  do i = 1,rstep

     dd1 = dd

     !***CALCULATE THE FREE FALL TIME***!
     rho = krome_get_rho(x(:))
     tff = sqrt(3d0 * krome_pi / (32d0*6.67e-8*rho))
     user_tff = tff
     dtH = 0.01d0 * tff        !TIME-STEP
     deldd = (dd/tff) * dtH
     dd = dd + deldd        !UPDATE DENSITY

     x(:) = x(:)*dd/dd1 !rescale abundances

     x(krome_idx_e) = krome_get_electrons(x(:))

     dt = dtH

     if(dd.gt.1d16) exit

     !local approximation for Av
     call krome_set_user_Av((dd*1d-3)**(2./3.))

     !solve the chemistry
     call krome(x(:),Tgas,dt)

     Tdust = krome_get_table_Tdust(x(:),Tgas)
     write(22,'(99E17.8e3)') dd, Tgas, Tdust, x(:)/dd

     !dump Tgas and normalized abundances
     !print every 100 steps
     if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas

  end do


  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
