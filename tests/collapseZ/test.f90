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
  integer::i,jz,imax(5)
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot,rho,zs(5)

  real*8::results(krome_nmols+3,10000,5)

  !INITIALIZE KROME PARAMETERS AND DUST
  call krome_init()
  call krome_set_zredshift(15d0)

  zs = (/-99d0, -4d0, -3d0, -2d0, -1d0/) !list of metallicities
  !$omp parallel do schedule(dynamic,1) default(none) &
  !$omp   private(jz,ntot,Tgas,x,dd,i,dd1,rho,tff,dt,dtH,deldd) &
  !$omp   shared(zs,imax,results)
  do jz = 1,size(zs)

     !INITIAL CONDITIONS
     ntot           = 0.1d0  !total density in 1/cm3
     Tgas           = 3d2    !temperature in kelvin

     !species initialization in 1/cm3
     x(:) = 1d-40

     x(KROME_idx_H)  = ntot          !H
     x(KROME_idx_H2) = 1d-6*ntot    !H2
     x(KROME_idx_E)  = 1d-4*ntot    !E
     x(KROME_idx_Hj) = 1d-4*ntot    !H+
     x(KROME_idx_HE) = 0.0775d0*ntot !He

     !rescale metallicity for neutral metals (C,Fe,Si,O)
     call krome_scale_Z(x(:), zs(jz))

     x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
     x(krome_idx_C)  = 1d-40

     !list abundances
     call krome_get_info(x(:),Tgas)

     dd = ntot

     print *,"solving..."
     print '(2a5,2a11)',"z","step","n(cm-3)","Tgas(K)"

     !loop over the hydro time-step
     do i = 1,rstep

        dd1 = dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3d0 * 3.1415d0 / (32d0*6.67d-8*rho))
        user_tff = tff
        dtH = 0.01d0 * tff !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd !UPDATE DENSITY

        x(:) = x(:) * dd / dd1 !rescale abundances

        dt = dtH

        if(dd>1d18) exit !quit after 1e18 1/cm3

        !solve the chemistry
        call krome(x(:),Tgas,dt)

        !dump Tgas and normalized abundances
        results(:,i,jz) = (/ zs(jz), dd, Tgas, x(:)/dd /)
        if(mod(i,100)==0) then
           print '(2I5,99E11.3)',jz,i,dd,Tgas !print every 100 steps
        end if
     end do
     imax(jz) = i - 1
  end do

  !dump all the results stored during the runs
  write(22,'(a)') "#Z rho Tgas "//krome_get_names_header()
  do jz = 1,size(zs)
     do i = 1,imax(jz)
        write(22,'(99E17.8e3)') results(:,i,jz)
     end do
     write(22,*)
     write(55,*)
  end do

  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
