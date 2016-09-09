!#######################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a primordial cloud.
!The dynamics is described by the Larson-Penston-type
! similar solution and includes cooling, heating processes
! and a UV background. For additional details look also to Omukai 2005.
!#######################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons

  integer,parameter::rstep = 500000
  integer::i,j,imax(5)
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot,rho,j21s(4)

  real*8::result(5,10000,4)

  !preset J21 values
  j21s = (/0d0, 1d0, 1d2, 1d5/)

  !INITIAL CONDITIONS
  call krome_set_zredshift(15d0) !redshift

  !INITIALIZE KROME PARAMETERS AND DUST
  call krome_init()

  !$omp parallel do schedule(dynamic,1) default(none) &
  !$omp   private(ntot,Tgas,x,dd,i,j,dd1,rho,tff,dt,dtH,deldd) &
  !$omp   shared(j21s,imax,result)
  do j=1,size(j21s)
     ntot           = 1d0    !total density in 1/cm3
     Tgas           = 1d2     !temperature in kelvin

     call krome_set_user_J21(j21s(j)) !common for J21

     !species initialization in 1/cm3
     x(:) = 1d-40

     x(KROME_idx_H)         = ntot           !H
     x(KROME_idx_H2)        = 1d-6*ntot    !H2
     x(KROME_idx_E)         = 1d-4*ntot    !E
     x(KROME_idx_Hj)        = 1d-4*ntot    !H+
     x(KROME_idx_HE)        = 0.0775*ntot    !He

     dd = ntot

     print *,"solving for J21=",krome_j21
     print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

     !loop over the hydro time-step
     do i = 1,rstep

        dd1 = dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3d0 * 3.1415d0 / (32d0*6.67d-8*rho))
        user_tff = tff
        dtH = 0.01d0 * tff        !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd        !UPDATE DENSITY

        x(:) = x(:) * dd / dd1

        dt = dtH

        if(dd.gt.1d16) exit

        !solve the chemistry
        call krome(x(:),Tgas,dt)

        !store results to dump later
        result(:,i,j) = (/ j21s(j),dd,Tgas,x(KROME_idx_H2)/dd,x(KROME_idx_H)/dd /)
        if(mod(i,100)==0) print '(2I5,99E11.3)',j,i,dd,Tgas

     end do
     imax(j) = i - 1
     print *,""
  end do

  !dump results
  write(22,'(a)') "#J21 ntot Tgas H2 H"
  do j = 1,size(j21s)
     do i = 1,imax(j)
        write(22,'(99E17.8e3)') result(:,i,j)
     end do
     write(22,*)
  end do

  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
