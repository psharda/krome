#include "fortran.def"
#include "phys_const.def"
#include "error.def"

  
  subroutine evaluate_temp(d, e, ge, u, v, w, &
    krome_x,imethod,idual,idim,tgas,temstart,utem,rhogas)
    use krome_user

    !     written by: KROME DEVELOPERS
    !     date: 2013
    !     
    !     PURPOSE:
    !     Evaluate temperature Tgas from the internal energy e
    !     or ge depending on the energy flag idual

    implicit none
#include "fortran_types.def"
    real*8::krome_x(:),p2d
    real*8::d,ge,e,tgas,rhogas
    real*8::u,v,w,temstart,utem,ntot
    integer::imethod,idim,idual,kgamma

    kgamma = krome_get_gamma(krome_x(:), Tgas)

    !Compute Pressure with various methods 
    if (imethod .eq. 2) then
       !Zeus - e() is really gas energy
       p2d = (kgamma - 1.d0)*d*e
    else
       if(idual.eq.1) then
          !PPM with dual energy -- use gas energy
          p2d = (kgamma - 1.d0)*d*ge
       else
          !PPM without dual energy -- use total energy
          p2d = e - 0.5d0*u**2
          if (idim.gt.1) p2d = p2d - 0.5d0 * v**2
          if (idim.gt.2) p2d = p2d - 0.5d0 * w**2
          p2d = max((kgamma - 1.d0) * d * p2d, tiny)
       endif
    endif

    !compute temperature
    tgas = max(p2d * utem / rhogas, temstart) 

  end subroutine evaluate_temp
