module krome_user_commons
  implicit none

#KROME_header
  
  !user can add here the common variables needed
  !for rate calculation (e.g. optical depth, CR rate, 
  !pressure, density, ...)
  real*8::tau,zrate,pah_size,gas_dust_ratio
  real*8::krome_J21,user_tff
  real*8::krome_invdvdz !inverse of |velocity gradient| 1/abs(dv/dz)

  !$omp threadprivate(krome_J21,user_tff,krome_invdvdz)

contains 

  !user can add here the functions he/she needs for
  !rate calculations (Kooij funcion provided as example)
  function kooij(kalpha,kbeta,kgamma,Tgas)
    real*8::kooij,kalpha,kbeta,kgamma,Tgas
    kooij = kalpha*(Tgas/3d2)**kbeta*exp(-kgamma/Tgas)
  end function kooij


end module krome_user_commons
