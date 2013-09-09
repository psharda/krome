!THIS IS THE ONE-ZONE COLLAPSE TEST
program test_krome

  use krome_commons
  use krome_main
  use krome_user
  use krome_user_commons


  integer,parameter::rstep = 500000
  integer::i,j
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(nmols),Tgas,dt
  real*8::ntot,rho,j21s(4)


  j21s = (/0.d0, 1d0, 1d2, 1d5/)
  do j=1,size(j21s)

     !INITIAL CONDITIONS
     redshift = 15d0    !redshift
     ntot = 1.d0           !total density in 1/cm3
     Tgas = 1d2              !temperature in kelvin

     !INITIALIZE KROME PARAMETERS AND DUST 
     call krome_init()

     krome_J21 = j21s(j)

     !species initialization in 1/cm3
     x(:) = 1.d-40

     x(KROME_idx_H)         = ntot           !H
     x(KROME_idx_H2)        = 1.0e-6*ntot    !H2
     x(KROME_idx_E)         = 1.0e-4*ntot    !E
     x(KROME_idx_Hj)        = 1.0e-4*ntot    !H+
     x(KROME_idx_HE)        = 0.0775*ntot    !He

     dd = ntot

     print *,"solving for J21=",krome_j21
     print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

     !loop over the hydro time-step
     do i = 1,rstep

        dd1 = dd

        !***CALCULATE THE FREE FALL TIME***!
        rho = krome_get_rho(x(:))
        tff = sqrt(3.0d0 * 3.1415d0 / (32.0d0*6.67e-8*rho))
        user_tff = tff
        dtH = 0.01d0 * tff        !TIME-STEP
        deldd = (dd/tff) * dtH
        dd = dd + deldd        !UPDATE DENSITY

        x(:) = x(:) * dd / dd1  

        dt = dtH 

        if(dd.gt.1d18) exit

        !solve the chemistry
        call krome(x(:),Tgas,dt)

        write(22,'(99E12.3e3)') j21s(j),dd,Tgas,x(KROME_idx_H2)/dd,x(KROME_idx_H)/dd
        if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas

     end do
     write(22,*)
     print *,""
  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
