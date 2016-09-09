!#########################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a primordial cloud irradiated
! by both a UV flux below the Lyman limit (13.6 eV) and a X-ray flux
! between 2-10 keV.
!For further details we refer to Inayoshi, MNRAS, 416, 2748-2759, 2011.
!Note that the results are different compared to Inayoshi+2011 as
! the cross-sections employed here for H and He photoionization are
! from Verner+1996, while Inayoshi adopted older outdated cross-sections.
!The dynamics is described by the Larson-Penston-type similar solution.
! For a direct benchmark see Fig. 3 of Latif+2014.
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
  real*8::ntot,rho,j21s(4)
  real*8::cool(12), heat(8),j21xs

  j21xs = 1d-1
  call krome_init()
  !preset J21 value for X-rays
  call krome_set_J21xray(j21xs)
  print *,"J21 Xray = ",j21xs
  !preset J21 values for UV flux
  j21s = (/0d0, 1d1, 1d3, 3d3/)
  do j=1,size(j21s)

     !INITIAL CONDITIONS
     call krome_set_zredshift(15d0)
     ntot = 1d-1            !total density in 1/cm3
     Tgas = 1.6d2            !temperature in kelvin

     !INITIALIZE KROME J21 parameter
     call krome_set_user_J21(j21s(j))

     !species initialization in 1/cm3
     x(:) = 1d-40

     x(KROME_idx_H)         = 0.9999*ntot  !H
     x(KROME_idx_H2)        = 2d-6*ntot    !H2
     x(KROME_idx_Hj)        = 2d-3*ntot    !H+
     x(KROME_idx_He)        = 0.0775*ntot  !He
     x(KROME_idx_E) = krome_get_electrons(x(:))

     dd = ntot

     print *,"solving for J21=",j21s(j)
     print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

     !output file header
     write(22,*) "#j21 ntot Tgas"

     !loop over the hydro time-step
     do i=1,rstep

        dd1 = dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3d0 * pi / (32d0*gravity*rho))
        user_tff = tff
        dtH = 0.01d0 * tff        !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd        !UPDATE DENSITY

        !rescale species number density to current density
        x(:) = x(:) * dd / dd1

        dt = dtH

        !run until 1e8 number density
        if(dd.gt.1d8) exit

        !solve the chemistry
        call krome(x(:),Tgas,dt)

        write(22,'(99E17.8e3)') j21s(j),dd,Tgas
        if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas

     end do
     write(22,*)

  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
