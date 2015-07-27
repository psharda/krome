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
  integer::i,jz
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot,rho,zs(5)

  !INITIALIZE KROME
  call krome_init()
  
  !black body (BB) temperatures (first is a dummy)
  zs = (/2.73d0, 2.73d0, 3d1, 3d2, 3d3/)
  do jz = 1,size(zs)
     print *,""
     print *,"****************"
     !INITIAL CONDITIONS
     call krome_set_zredshift(0d0)
     ntot           = 0.1d0  !total density in 1/cm3
     Tgas           = 3d2    !temperature in kelvin

     !species initialization in 1/cm3
     x(:) = 1.d-40

     x(KROME_idx_Hj)  = ntot          !H
     x(KROME_idx_H2)  = 1.d-6*ntot    !H2
     x(KROME_idx_HEj) = 0.0775d0*ntot !He

     !rescale metallicity for neutral metals (C,Fe,Si,O)
     call krome_scale_Z(x(:), -3d0)

     x(krome_idx_Cj)  = x(krome_idx_C)  !carbon is fully ionized
     x(krome_idx_C)   = 1d-40
     x(krome_idx_Sij) = x(krome_idx_Si) !silicon is fully ionized
     x(krome_idx_Si)  = 1d-40
     x(krome_idx_Fej) = x(krome_idx_Fe) !iron is fully ionized
     x(krome_idx_Fe)  = 1d-40

     x(krome_idx_E) = krome_get_electrons(x)

     !init BB radiation from temperature with automatic limits
     call krome_set_photoBin_BBlog_auto(zs(jz))

     !first is primordial
     if(jz==1) call krome_photoBin_scale(0d0)
     dd = ntot

     print *,"solving..."
     print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

     !loop over the hydro time-step
     do i = 1,rstep

        dd1=dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3.0d0 * 3.1415d0 / (32.0d0*6.67e-8*rho))
        user_tff = tff
        dtH = 0.01d0 * tff        !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd        !UPDATE DENSITY

        x(:) = x(:)*dd/dd1 !rescale abundances

        dt = dtH 

        if(dd.gt.1d12) exit !quit after 1e16 1/cm3

        x(krome_idx_e) = krome_get_electrons(x)
        !solve the chemistry
        call krome(x(:),Tgas,dt)

        !dump Tgas and normalized abundances
        write(22,'(I8,99E17.8e3)') jz,zs(jz),dd,Tgas,x(:)/dd 
        if(mod(i,100)==0) then 
           print '(I5,99E11.3)',i,dd,Tgas !print every 100 steps
        end if
        write(88,'(999E17.8e3)') dd,Tgas,krome_get_cooling_array(x,Tgas)
        write(77,'(999E17.8e3)') dd,Tgas,krome_get_heating_array(x,Tgas)
        write(66,'(999E17.8e3)') dd,Tgas,x(:)
        !if(dd>1e8) exit
     end do
     write(22,*)
     write(66,*)
     write(77,*)
     write(88,*)
  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
