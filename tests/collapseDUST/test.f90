!THIS IS THE ONE-ZONE COLLAPSE TEST
program test_krome

  use krome_commons
  use krome_main
  use krome_user
  use krome_user_commons
  use krome_dust


  integer,parameter::rstep = 500000,nd=1,ndtype=1
  integer::i,j
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(nmols),Tgas,dt
  real*8::ntot,rho, dust_gas_ratio(ndtype)
  real*8::xdust(nd),adust(nd),js(4)

  js(:) = (/0.d0, 1d-6, 1d-4, 1d-2/)
  do j=1,size(js)
     !INITIAL CONDITIONS
     redshift = 0d0    !redshift
     ntot = 1.d0           !total density in 1/cm3
     Tgas = 3d2              !temperature in kelvin

     !species initialization in 1/cm3
     x(:) = 1.d-40

     x(KROME_idx_H)         = ntot           !H
     x(KROME_idx_H2)        = 1.0e-6*ntot    !H2
     x(KROME_idx_E)         = 1.0e-4*ntot    !E
     x(KROME_idx_Hj)        = 1.0e-4*ntot    !H+
     x(KROME_idx_HE)        = 0.0775*ntot    !He

     dust_gas_ratio(:) = js(j)

     !INITIALIZE KROME PARAMETERS AND DUST 
     call krome_init()
     call krome_init_dust(xdust(:),adust(:),dust_gas_ratio(:),x(:))

     dd = ntot

     print *,xdust(:),dust_gas_ratio(:)
     print *,"solving..."
     print '(a5,3a11)',"step","n(cm-3)","Tgas(K)","Tdust(K)"

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

        x(:) = x(:) * dd/dd1  
        xdust(:) = xdust(:) * dd/dd1
        dt = dtH 

        if(dd.gt.1d18) exit

        !solve the chemistry
        call krome(x(:),Tgas,dt,xdust(:))

        write(22,'(I5,99E12.3e3)') j,js(j),dd,Tgas,krome_dust_T(:),&
             x(:)/dd,xdust(:)/dd

        if(mod(i,100)==0) print '(I5,99E11.3)',i,dd,Tgas,krome_dust_T(:)
        
     end do
     write(22,*)
  end do
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
