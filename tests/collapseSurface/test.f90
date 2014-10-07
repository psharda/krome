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
  integer,parameter::nz=3
  integer,parameter::rstep = 500000
  integer::i,jz
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot,rho,zs(nz),zred

  zred = 0d0
  !INITIALIZE KROME PARAMETERS AND DUST 
  call krome_init()
  call krome_set_zredshift(zred)
  call krome_set_Tcmb(2.73d0*(zred+1d0))
  zs = (/-5d0, -4d0, -2d0/) !list of metallicities

  do jz = 1,size(zs)

     call krome_set_dust_distribution()

     !INITIAL CONDITIONS
     ntot           = 0.1d0  !total density in 1/cm3
     Tgas           = 3d2    !temperature in kelvin

     !species initialization in 1/cm3
     x(:) = 1.d-40

     x(KROME_idx_H)  = ntot          !H
     x(KROME_idx_H2) = 1.d-6*ntot    !H2
     x(KROME_idx_E)  = 1.d-4*ntot    !E
     x(KROME_idx_Hj) = 1.d-4*ntot    !H+
     x(KROME_idx_HE) = 0.0775d0*ntot !He

     call krome_scale_dust_gas_ratio(1d-2*1d1**zs(jz),x(:))
     call krome_set_defaultTdust((zred+1d0)*2.73d0)

     !rescale metallicity for neutral metals (C,Fe,Si,O)
     call krome_scale_Z(x(:), zs(jz))

     x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
     x(krome_idx_C)  = 1d-40

     !list abundances
     !call krome_get_info(x(:),Tgas)

     dd = ntot

     print *,"solving..."
     print '(2a5,2a11)',"z","step","n(cm-3)","Tgas(K)","Tdust(K)"

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

        x(krome_idx_e) = krome_get_electrons(x(:))
        dt = dtH 

        if(dd>1d18) exit !quit after 1e10 1/cm3
        call krome_scale_dust_distribution(dd/dd1)

        !dust evaporation
        if(Tgas>1.5d3) call krome_scale_dust_distribution(0d0)

        !solve the chemistry
        call krome(x(:),Tgas,dt)

        !dump Tgas and normalized abundances
        write(66,'(I5,99E17.8e3)') jz,dd,Tgas,krome_get_Tdust()
        if(mod(i,100)==0) then
           print '(2I5,99E11.3)',jz,i,dd,Tgas,krome_get_Tdust()
        end if
     end do
     write(33,*)
     write(34,*)
     write(66,*)
  end do


  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
