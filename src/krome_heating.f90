module KROME_heating
contains

#KROME_header

  !*******************************
  function heating(n, Tgas)
    implicit none
    real*8::n(:), Tgas
    real*8::heating 
    !total heating erg/cm3/s
    heating = 0.d0
#IFKROME_useHeatingChem
    heating = heating + heatingChem(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingCompress
    heating = heating + heat_compress(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhoto
    heating = heating + photo_heating(n(:))
#ENDIFKROME

#IFKROME_useHeatingdH
    heating = heating + heat_dH(n(:),Tgas)
#ENDIFKROME
  end function heating
  

#IFKROME_useHeatingdH
  !*************************
  function heat_dH(n,Tgas)
    !heating from reaction enthalpy erg/s/cm3
    use krome_commons
    implicit none
    real*8::heat_dH,heat,n(:),Tgas,T4
    real*8::logT,lnT,Te,lnTe,T32,t3,invT,invTe,sqrTgas,invsqrT32,sqrT32
#KROME_vars
    
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
    n0 = 1d99 !density for fake opacity (see KROME paper)
    photo_heating = 0.d0
#KROME_photo_heating

    photo_heating = photo_heating * eV_to_erg
  end function photo_heating
#ENDIFKROME

#IFKROME_useHeatingChem
  !H2 FORMATION HEATING
  !UNITS = erg/cm3/s
  !Following Hollenbach & McKee (1979), we assume the heat deposited 
  !per a formed H2 is: 
  !3.53*(1+ncr/n)^−1 eV for H2 formation by H− process
  !1.83*(1+ncr/n)^−1 eV for H2 formation by H2+ process
  !4.48*(1+ncr/n)^−1 eV for H2 formation by the three-body 
    !*******************************
  function heatingChem(n, Tgas)
    use krome_constants
    use krome_commons
    implicit none
    real*8::heatingChem, n(:), Tgas
    real*8::h2heatfac,HChem,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2,dd
    real*8::logT,lnT,Te,lnTe,T32,invT,invTe,sqrTgas,invsqrT32,sqrT32
    real*8::Tgas2,Tgas3,Tgas4,T0,T02,T03,T04,T0inv, t
    dd = sum(n(1:nmols))

#KROME_Tshortcuts

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**(2.d0))
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    
#KROME_HChem_terms

    HChem = HChem * h2heatfac
    
    heatingChem = HChem * eV_to_erg  !erg/cm3/s
    
  end function heatingChem
#ENDIFKROME

#IFKROME_useHeatingCompress
  !***********************
  function heat_compress(n, Tgas)
    use krome_user_commons
    use krome_commons
    use krome_constants
    real*8::heat_compress,n(:), dd, Tgas

    dd = sum(n(1:nmols)) !total number density

    !COMPRESSIONAL HEATING
    !note that krome_gamma is computed in the FEX in krome_ode module
    heat_compress = dd * boltzmann_erg * Tgas * (krome_gamma - 1.d0) / user_tff !erg/s/cm3

  end function heat_compress
#ENDIFKROME
  
end module KROME_heating
