!################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a primordial cloud.
!The dynamics is described by the Larson-Penston-type
! similar solution and includes cooling and heating processes.
!For additional details look also to Omukai 2000 and
! the KROME paper.
!################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons


  integer,parameter::rstep = 500000
  integer::i
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot,rho

  !INITIAL CONDITIONS
  krome_redshift = 15d0    !redshift
  ntot = 1d0                !total density in 1/cm3
  Tgas = 1d2                 !temperature in kelvin

  !INITIALIZE KROME PARAMETERS AND DUST
  call krome_init()

  !species initialization in 1/cm3
  x(:) = 1d-40

  x(KROME_idx_H)         = ntot           !H
  x(KROME_idx_H2)        = 1d-6*ntot    !H2
  x(KROME_idx_E)         = 1d-4*ntot    !E
  x(KROME_idx_Hj)        = 1d-4*ntot    !H+
  x(KROME_idx_HE)        = 0.0775*ntot    !He

  dd = ntot

  print *,"solving..."
  print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

  !output header
  write(22,*) "#ntot Tgas "//trim(krome_get_names_header())

  !loop over the hydro time-step
  do i = 1,rstep

     dd1=dd

     !***CALCULATE THE FREE FALL TIME***!
     rho = krome_get_rho(x(:))
     tff = sqrt(3d0 * 3.1415d0 / (32d0*6.67d-8*rho))
     user_tff = tff
     dtH = 0.01d0 * tff        !TIME-STEP
     deldd = (dd/tff) * dtH
     dd = dd + deldd        !UPDATE DENSITY

     x(:) = x(:)*dd/dd1

     dt = dtH

     if(dd.gt.1d18) exit

     !solve the chemistry
     call krome(x(:),Tgas,dt)

     write(22,'(99E17.8e3)') dd,Tgas,x(:)/dd
     if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas

  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
