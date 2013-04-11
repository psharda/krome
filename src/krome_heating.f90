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
#IFKROME_useHeatingA
    heating = heating + heatingA(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingCompress
    heating = heating + heat_compress(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhoto
    heating = heating + photo_heating(n(:))
#ENDIFKROME
    
  end function heating
  
  
#IFKROME_useHeatingPhoto
  !**************************
  function photo_heating(n)
    use krome_commons
    real*8::photo_heating,n(:),n0
    n0 = 1d-1 !density for fake opacity (see KROME paper)
    photo_heating = 0.d0
#KROME_photo_heating
  end function photo_heating
#ENDIFKROME

#IFKROME_useHeatingA
  !H2 FORMATION HEATING
  !OMUKAI 2000
  !UNITS=erg/cm3/s
  !Following Hollenbach & McKee (1979), we assume the heat deposited 
  !per a formed H2 is: 
  !3.53*(1+ncr/n)^−1 eV for H2 formation by H− process (k_HM) 
  !1.83*(1+ncr/n)^−1 eV for H2 formation by H2+ process (k_H2j)
  !4.48*(1+ncr/n)^−1 eV for H2 formation by the three-body 
  !reactions (k_3B, k_3B1)
  !NOTE: THE LAST TERM IS A COOLING SOURCE TERM FROM k_3B2 COSA VUOI FARE?
  !*******************************
  function heatingA(n, Tgas)
    use krome_constants
    use krome_commons
    implicit none
    real*8::heatingA, n(:), Tgas
    real*8::k_3B2, k_HM, k_H2j, k_3B, k_3B1
    real*8::h2heatfac,H2delta,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2
    real*8::Te,invTe,dd

    dd = sum(n(1:nmols))

    Te = Tgas*boltzmann_eV
    invTe = 1.d0/te

    heatingA = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**(2.d0))
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    k_3B  = 5.5d-29/Tgas  !3H --> H2 + H
    k_3B1 = k_3B/8.0d0    !2H + H2 --> H2 + H 
    k_3B2 = (6.5d-7/Tgas**0.5d0)*exp(-5.2d4/Tgas)&
         *(1.d0-exp(-6.0d3/Tgas)) !H2 + H --> 3H

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    H2delta = n(idx_H)*(4.48d0*k_3B*n(idx_H)**(2.d0) &
         + 4.48d0*(k_3B1)*n(idx_H)*n(idx_H2)*2.0 &
                                !+ 3.53d0*k_HM*n(idx_Hk) &
                                !+ 1.83*k_H2j*n(idx_H2)
         )*h2heatfac &
         - n(idx_H)*(4.48d0*k_3B2*n(idx_H2))


    if(H2delta.gt.0.0d0)then
       heatingA = H2delta*eV_to_erg  !erg/cm3/s
    else
       heatingA = 0.0d0
    endif

  end function heatingA
#ENDIFKROME

#IFKROME_useHeatingCompress
  !***********************
  function heat_compress(n, Tgas)
    use krome_user_commons
    use krome_commons
    use krome_constants
    real*8::heat_compress,n(:), eint, gamma
    real*8::dd, p2d, Tgas
    gamma = 5./3.
    dd = sum(n(1:nmols))
    eint = (1.0/(gamma-1.0d0))*(boltzmann_erg*Tgas)/(1.22*p_mass) !erg/g 
    p2d = (gamma-1.0d0)*dd*eint  !erg/g/cm3

    !COMPRESSIONAL HEATING
    heat_compress = (p2d)*(p_mass/user_tff) !erg/s/cm3
  end function heat_compress
#ENDIFKROME
  
end module KROME_heating
