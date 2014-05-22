module KROME_heating
contains

#KROME_header

  !************************
  function heating(n,Tgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::heating
    
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))
    
  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(8),heats(8)
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

#IFKROME_useHeatingPhotoAv
    heats(5) = heat_photoAv(n(:),Tgas,k(:))
#ENDIFKROME

#IFKROME_useHeatingCR
    heats(6) = heat_CR(n(:),Tgas,k(:))
#ENDIFKROME

#IFKROME_useHeatingPhotoDust
    heats(7) = heat_photoDust(n(:),Tgas)
#ENDIFKROME

#IFKROME_useHeatingXRay
    heats(8) = heat_XRay(n(:),Tgas,k(:))
#ENDIFKROME
    
    get_heating_array(:) = heats(:)

    !remove the comment below to write heating terms to fort.55
    !write(55,'(99E17.8e3)') sum(n(1:nmols)),Tgas,heats(:)

    !gnuplot command (n=100, and m=1 for density or m=2 for temperature) 
    !plot 'fort.55' u m:(abs($3)) every n w l t "chem",\
    ! '' u m:4 every n w l t "compress",\
    ! '' u m:5 every n w l t "photo",\
    ! '' u m:6 every n w l t "enthalpy"

  end function get_heating_array

#IFKROME_useHeatingXRay
  !*************************
  !heating from xrays in erg/s/cm3
  function heat_XRay(n,Tgas,k)
    use krome_commons
    use krome_subs
    implicit none
    real*8::n(:),Tgas,heat_Xray,k(:),ntot
    real*8::xheat_H,xheat_He,logH,logHe
    real*8::xe,ratexH,ratexHe

    ntot = get_Hnuclei(n(:))
    xe = min(n(idx_e)/ntot,1d0)

    logH = log10(max(n(idx_H),1d-40))
    logHe = log10(max(n(idx_He),1d-40))

    heat_Xray = 0d0
    xheat_H = fit_anytab2D(user_xheat_H_anytabx(:), &
         user_xheat_H_anytaby(:), &
         user_xheat_H_anytabz(:,:), &
         user_xheat_H_anytabxmul, &
         user_xheat_H_anytabxmul, &
         logH,logHe-logH)
    xheat_He = fit_anytab2D(user_xheat_He_anytabx(:), &
         user_xheat_He_anytaby(:), &
         user_xheat_He_anytabz(:,:), &
         user_xheat_He_anytabxmul, &
         user_xheat_He_anytabxmul, &
         logH,logHe-logH)

    !prepares varibles for xray photochemistry
    ratexH = 1d1**xheat_H
    ratexHe = 1d1**xheat_He
    
    heat_Xray = ratexH * n(idx_H)
    heat_Xray = ratexHe * n(idx_He)
    heat_Xray = heat_Xray * .9971d0 * (1d0-(1d0*xe**.2663)**1.3163)

  end function heat_XRay
#ENDIFKROME

#IFKROME_useHeatingPhotoDust
  !***************************
  function heat_photoDust(n,Tgas)
    !photoelectric effect from dust in erg/s/cm3
    use krome_commons
    use krome_subs
    implicit none
    real*8::heat_photoDust,n(:),Tgas,ntot,eps
    real*8::Ghab,z,zsun,psi

    ntot = get_Hnuclei(n(:))
    zsun = 0.02d0 !solar metallicity
    Ghab = 1.69d0 !habing flux, 1.69 is Draine78
    psi = 2.d0*Ghab*sqrt(Tgas)*n(idx_e)
    eps = 4.9d-2/(1d0+4d-3*psi**.73) + &
         3.7d-2*(Tgas*1d-4)**.7/(1d0+2d-4*psi)
    z = #KROME_photoDustZ !metallicty
    heat_photoDust = 1.3d-24*eps*Ghab*ntot*z/zsun

  end function heat_photoDust
#ENDIFKROME

#IFKROME_useHeatingPhotoAv
  !******************************
  function heat_photoAv(n,Tgas,k)
    !heating from  photoreactions using rate approximation erg/s/cm3
    use krome_commons
    use krome_user_commons
    use krome_subs
    implicit none
    real*8::heat_photoAv,n(:),Tgas,k(:)
    real*8::ncrn,ncrd1,ncrd2,yH,yH2,ncr,h2heatfac,dd,Rdiss

    dd = get_Hnuclei(n(:))
    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))
    
    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless
    
    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless
    
    Rdiss = #KROME_RdissH2

    !photodissociation H2 heating
    heat_photoAv = 6.4d-13*Rdiss*n(idx_H2)

    !UV photo-pumping H2
    heat_photoAv = heat_photoAv + 2.7d-11*Rdiss*h2heatfac*n(idx_H2)
    
  end function heat_photoAv
#ENDIFKROME

#IFKROME_useHeatingCR
  !***************************
  function heat_CR(n,Tgas,k)
    !heating from cosmic rays erg/s/cm3
    use krome_commons
    implicit none
    real*8::heat_CR,n(:),Tgas,Hfact,k(:)

    Hfact = 3.20435313d-11 !erg

    heat_CR = 0d0

#KROME_heatingCR

  end function heat_CR
#ENDIFKROME

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
    !photo heating in erg/cm3/s using bin-based
    ! approach. Terms are computed in the
    ! krome_photo module
    use krome_commons
    use krome_constants
    implicit none
    real*8::photo_heating,n(:)

    photo_heating = 0.d0
#KROME_photo_heating

  end function photo_heating
#ENDIFKROME

#IFKROME_useHeatingChem
  !H2 FORMATION HEATING and other exo/endothermic 
  ! processes (including H2 on dust) in erg/cm3/s
  !krome builds the heating/cooling term according
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
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
    dd = get_Hnuclei(n(:))

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = #KROME_small

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
