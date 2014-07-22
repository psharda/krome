!#########################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a primordial cloud irradiated
! by both a UV flux below the Lyman limit (13.6 eV) and a X-ray flux
! between 2-10 keV. 
!For further details we remind to Inayoshi, MNRAS, 416, 2748-2759, 2011.
!Note that the results are different compared to Inayoshi+, 2011 as
! the cross-sections employed here for H and He photoionization are
! from Verner+ 1996, while Inayoshi+ adopted older outdated cross-sections.
!The dynamics is described by the Larson-Penston-type similar solution.
!##########################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_constants
  use krome_subs
  use krome_photo

  integer,parameter::rstep = 500000
  integer::i,j
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt,f(krome_nmols)
  real*8::ntot,rho,j21s(7),mass(krome_nspec)
  real*8::cool(12), heat(8)
  

  !preset J21 value for X-rays
  user_krome_J21xr = 1.0d-2
  !preset J21 values for UV flux
  j21s = (/0.0d0, 1d-1, 5d0, 1d1, 2d1, 5d1, 1d2/)
  do j=1,size(j21s)

     !INITIAL CONDITIONS
     krome_redshift = 15d0   !redshift
     ntot = 1.d-1            !total density in 1/cm3
     Tgas = 1.6d2            !temperature in kelvin

     !INITIALIZE KROME PARAMETERS AND DUST 
     krome_J21 = j21s(j) !common for J21
     call krome_init()

     !species initialization in 1/cm3
     x(:) = 1.d-40

     x(KROME_idx_H)         = 0.9999*ntot    !H
     x(KROME_idx_H2)        = 2.0e-6*ntot    !H2
     x(KROME_idx_E)         = 2.0e-3*ntot    !1.0e-4*ntot    !E
     x(KROME_idx_Hj)        = 2.0e-3*ntot    !1.0e-4*ntot    !H+
     x(KROME_idx_HE)        = 0.0775*ntot    !He


     mass(:) = get_mass()
     dd = ntot

     print *,"solving for J21=",krome_j21
     print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

     !loop over the hydro time-step
     do i = 1,rstep

        dd1 = dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3.0d0 * pi / (32.0d0*gravity*rho))
        user_tff = tff
        dtH = 0.01d0 * tff        !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd        !UPDATE DENSITY

        x(:) = x(:) * dd / dd1  

        dt = dtH 

        if(dd.gt.1d8) exit

        !store the heating/cooling contributions
        cool = krome_get_cooling_array(x, Tgas)
        heat = krome_get_heating_array(x, Tgas)


        !solve the chemistry
        call krome(x(:),Tgas,dt)

        write(22,'(99E17.8e3)') j21s(j),dd,Tgas
        write(23,'(99E17.8e3)') j21s(j),dd, cool
        write(24,'(99E17.8e3)') j21s(j),dd, heat 
        if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas

     end do
     write(22,*)
     print *,""
  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
