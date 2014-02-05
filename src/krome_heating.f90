module KROME_heating
contains

#KROME_header

  !*******************************
  function heating(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::heating,heats(4)
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

#IFKROME_useHeatingChem
    heats(1) = heatingChem(n(:), Tgas, k(:), nH2dust)
#ENDIFKROME

#IFKROME_useHeatingCompress
    heats(2) = heat_compress(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhoto
    heats(3) = photo_heating(n(:))
#ENDIFKROME

#IFKROME_useHeatingdH
    heats(4) = heat_dH(n(:),Tgas)
#ENDIFKROME
    
    heating = sum(heats)

    !remove the comment below to write heating terms to fort.55
    !write(55,'(99E17.8e3)') sum(n(1:nmols)),Tgas,heats(:)

    !gnuplot command (n=100, and m=1 for density or m=2 for temperature) 
    !plot 'fort.55' u m:(abs($3)) every n w l t "chem",\
    ! '' u m:4 every n w l t "compress",\
    ! '' u m:5 every n w l t "photo",\
    ! '' u m:6 every n w l t "enthalpy"

  end function heating
  

#IFKROME_useHeatingdH
  !*************************
  function heat_dH(n,Tgas)
    !heating from reaction enthalpy erg/s/cm3
    use krome_commons
    implicit none
    real*8::heat_dH,heat,n(:),Tgas,T4,small
    real*8::logT,lnT,Te,lnTe,T32,t3,invT,invTe,sqrTgas,invsqrT32,sqrT32
#KROME_vars
    small = 1d-40
    logT = log10(Tgas) !log10 of Tgas (#)
    lnT = log(Tgas) !ln of Tgas (#)
    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    T32 = Tgas/3.d2 !Tgas/(300 K) (#)
    t3 = T32 !alias for T32 (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    invTe = 1.d0/Te !inverse of T (1/eV)
    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)
    invsqrT32 = 1.d0/sqrt(T32)
    sqrT32 = sqrt(T32)

    heat = 0.d0

#KROME_rates
#KROME_dH_heating

    heat_dH = heat    

  end function heat_dH
#ENDIFKROME
  
#IFKROME_useHeatingPhoto
  !**************************
  function photo_heating(n)
    use krome_commons
    use krome_constants
    real*8::photo_heating,n(:),n0
    n0 = 1d99 !density for fake opacity
    photo_heating = 0.d0
#KROME_photo_heating

    photo_heating = photo_heating * eV_to_erg
  end function photo_heating
#ENDIFKROME

#IFKROME_useHeatingChem
  !H2 FORMATION HEATING
  !UNITS = erg/cm3/s
  !krome build the heating/cooling term according
  ! to the chemical network employed
  !*******************************
  function heatingChem(n, Tgas, k, nH2dust)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_subs
    implicit none
    real*8::heatingChem, n(:), Tgas,k(:),nH2dust
    real*8::h2heatfac,HChem,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small
    dd = get_Hnuclei(n(:))

    small = 1d-40

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    n2H = n(idx_H) * n(idx_H)

#KROME_HChem_terms
#KROME_HChem_dust

    heatingChem = HChem * eV_to_erg  !erg/cm3/s

  end function heatingChem
#ENDIFKROME

#IFKROME_useHeatingCompress
  !***********************
  !evaluates compressional heating
  ! WARNING: user_tff is a common variable
  ! available in krome_user_commons.f90
  function heat_compress(n, Tgas)
    use krome_user_commons
    use krome_commons
    use krome_constants
    use krome_subs
    real*8::heat_compress,n(:), dd, Tgas

    dd = get_Hnuclei(n(:)) !sum(n(1:nmols)) !total number density

    !COMPRESSIONAL HEATING
    heat_compress = dd * boltzmann_erg * Tgas / user_tff !erg/s/cm3

  end function heat_compress
#ENDIFKROME
  
end module KROME_heating
