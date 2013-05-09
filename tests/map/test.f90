program test
  use krome_main
  use krome_photo
  use krome_user
  use krome_user_commons
  implicit none
  integer,parameter::nsp=12
  real*8::x(nsp),rho,Tgas,dt,spy,Tgaso,dtmax
  real*8::Tmin,Tmax,rhomin,rhomax,xo(nsp)
  integer::irmax,itmax,it,ir
  

  call krome_init()
  call krome_init_photo()

  itmax = 10
  irmax = 10
  Tmin = log10(1.d2)
  Tmax = log10(1.d6)
  rhomin = log10(1.d-24)
  rhomax = log10(1.d-16)
  dtmax = 1d6

  xo(:) = 1d-20
  xo(KROME_idx_H)     = 0.9225d0  !H
  xo(KROME_idx_E)     = 1.0d-4    !E
  xo(KROME_idx_Hj)    = 1.0d-4    !H+
  xo(KROME_idx_D)     = 1.0d-20   !D
  xo(KROME_idx_Dj)    = 1.0d-20   !D+
  xo(KROME_idx_HE)    = 0.0972d0  !He
  xo(KROME_idx_HEj)   = 1.0d-20   !He+
  xo(KROME_idx_H2j)   = 1.0d-20   !H2+
  xo(KROME_idx_H2)    = 1.0d-5    !H2
  xo(KROME_idx_HD)    = 1.0d-8    !HD
  xo(KROME_idx_Hk)    = 1.0d-20   !H-
  xo(KROME_idx_HEjj)  = 1.0d-20   !He++
  xo(:) = xo(:) / sum(xo)
  spy = 365d0 * 24d0 * 3600d0
  print '(2a12)',"Tgas (K)","rho (g/cm3)"
  do it=itmax,1,-1
     do ir=irmax,1,-1
        x(:) = xo(:)
        dt = spy * dtmax
        rho = 1d1**((irmax-ir)*(rhomax-rhomin)/(irmax-1) + rhomin) !g/cm3
        Tgas = 1d1**((itmax-it)*(Tmax-Tmin)/(itmax-1) + Tmin) !g/cm3
        Tgaso = Tgas

        print '(2E12.3)',Tgaso,rho

        call krome(x(:), rho, Tgas, dt)
        write(66,'(99E12.3e3)') rho,Tgaso,Tgas,x(:)
     end do !end rho loop
     write(66,*)
  end do !end Tgas loop
  print *,"End of map test!"
  print *,"plot.gps in gnuplot to plot the map!"
end program test
