!##################################################
!This test represents a series of one-zone models
! with different photon flux (J21) and an initial
! temperature of 10^4 Kelvin.
!##################################################
program test
  use krome_main
  use krome_photo
  use krome_user
  use krome_user_commons
  use krome_constants
  implicit none
  real*8::x(krome_nmols),rho,Tgas,dt,Tgaso,dtmax
  real*8::J21min,J21max,rhomin,rhomax,xo(krome_nmols)
  integer::irmax,itmax,it,ir

  call krome_init()

  itmax = 50 !number of J21 values
  irmax = 50 !number of rho values
  J21min = log10(1.d-3) !min J21 value (log)
  J21max = log10(1.d1) !max J21 value (log)
  rhomin = log10(1.d-25) !min rho value (log)
  rhomax = log10(1.d-21) !max rho value (log)

  xo(:) = 1d-20 !default value
  xo(KROME_idx_H)     = 0.9225d0  !H
  xo(KROME_idx_E)     = 1.0d-4    !E
  xo(KROME_idx_Hj)    = 1.0d-4    !H+
  xo(KROME_idx_HE)    = 0.0972d0  !He
  xo(KROME_idx_H2)    = 1.0d-5    !H2
  xo(KROME_idx_HD)    = 1.0d-8    !HD

  Tgaso = 1d4 !initial temperature
  xo(:) = xo(:) / sum(xo) !normalize

  print '(2a12)',"J21"
  !loop on J21s
  do it=1,itmax
     krome_J21 = 1d1**((it-1)*(J21max-J21min)/(itmax-1) + J21min) !J21 
     call krome_init_photo()
     print '(E12.3,F12.2,a1)',krome_J21,it*1d2/itmax,"%" 
     !loop on densities
     do ir=1,irmax
        x(:) = xo(:) !copy initial species
        rho = 1d1**((ir-1)*(rhomax-rhomin)/(irmax-1) + rhomin) !g/cm3
        dt = sqrt(3.d0*pi/32.d0/gravity/rho) !free-fall time
        Tgas = Tgaso !copy initial Tgas
        call krome(x(:), rho, Tgas, dt) !KROME
        write(66,'(99E12.3e3)') rho,krome_J21,Tgaso,Tgas,x(:) !dump
     end do !end rho loop
     write(66,*)
  end do !end J21 loop
  print *,"Models done:",itmax*irmax
  print *,"End of map test!"
  print *,"plot.gps in gnuplot to plot the map!"

end program test
