module krome_user_commons
  implicit none

#KROME_header
  
  !user can add here the common variables needed
  !for rate calculation (e.g. optical depth, CR rate, 
  !pressure, density, ...)
  real*8::Pgas,rhogas
  real*8::tau,zrate,pah_size,gas_dust_ratio

contains 

  !user can add here the functions he/she needs for
  !rate calculations (Kooij funcion provided as example)
  function kooij(kalpha,kbeta,kgamma,Tgas)
    real*8::kooij,kalpha,kbeta,kgamma,Tgas
    kooij = kalpha*(Tgas/3d2)**kbeta*exp(-kgamma/Tgas)
  end function kooij

  !******************
  function ktrimol(zazero,zaone,zcfirst,zcsecond,Tgas,dn)
    real*8::ktrimol,zazero,zaone,zcfirst,zcsecond,Tgas,dn
    ktrimol = 0.d0
    !ktrimol = (zazero*(3d2/Tgas)**zcfirst) &
    !     *dn / (1.d0+(zazero*(3d2/Tgas)**zcfirst)*dn &
    !     /(zaone*(3d2/Tgas)**zcsecond))&
    !     *0.6d0**(1.d0/(1.d0 + (log10((zazero*(3d2/Tgas)**zcfirst)&
    !    *dn/(zaone*(3d2/Tgas)**zcsecond)))**2) )
  end function ktrimol
end module krome_user_commons
