module KROME_heating
contains

#KROME_header

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

    f2 = 1.
#IFKROME_useHeatingGH
    !this parameter controls the smoothness of the
    ! merge between the two cooling / heating functions
    smooth = 1.d-3

    !smoothing functions | f1+f2=1
    !f1 = (tanh(smooth*(Tgas-1d4))+1.d0)*0.5d0
    f2 = (tanh(smooth*(-Tgas+1d4))+1.d0)*0.5d0

    !heating is already included in the GH cooling function (that thus can be negative), and f1 is not needed
#ENDIFKROME

#IFKROME_useHeatingChem
    heats(idx_heat_chem) = heatingChem(n(:), Tgas, k(:), nH2dust)
#ENDIFKROME

#IFKROME_useHeatingCompress
    heats(idx_heat_compress) = heat_compress(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhoto
    heats(idx_heat_photo) = photo_heating(n(:),f2)
#ENDIFKROME

#IFKROME_useHeatingdH
    heats(idx_heat_dH) = heat_dH(n(:),Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhotoAv
    heats(idx_heat_photoAv) = heat_photoAv(n(:),Tgas,k(:))
#ENDIFKROME

#IFKROME_useHeatingCR
    heats(idx_heat_CR) = heat_CR(n(:),Tgas,k(:))
#ENDIFKROME

#IFKROME_useHeatingPhotoDust
    heats(idx_heat_dust) = heat_photoDust(n(:),Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhotoDustNet
    heats(idx_heat_dust) = heat_netPhotoDust(n(:),Tgas)
#ENDIFKROME

#IFKROME_useHeatingXRay
    heats(idx_heat_xray) = heat_XRay(n(:),Tgas,k(:))
#ENDIFKROME

#IFKROME_useHeatingVisc
    heats(idx_heat_visc) = heat_Visc(n(:),Tgas)
#ENDIFKROME


#IFKROME_useHeatingXRay
    heats(idx_heat_xray) = f2 * heats(idx_heat_xray)
#ENDIFKROME

#IFKROME_useHeatingZExtended
    !this parameter controls the smoothness of the
    ! merge between the two cooling / heating functions
    smooth = 1.d-3
    !smoothing functions | f1+f2=1 | f1 is not needed
    f2 = = (tanh(smooth*(Tgas-1d4))+1.d0)*0.5d0 
#ENDIFKROME

#IFKROME_useHeatingZCIE
    heats(idx_heat_ZCIE) = f2 * heat_ZCIE(n(:),Tgas)
#ENDIFKROME

    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array


  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom
#KROME_custom_heating_var_define

    heat_custom = 0d0
#KROME_custom_heating_var
#KROME_custom_heating_expr

  end function heat_custom


#IFKROME_useHeatingZCIE
  function heat_ZCIE(n,inTgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer,parameter::imax=coolZCIEn1
    integer,parameter::jmax=coolZCIEn2
    integer,parameter::kmax=coolZCIEn3
    integer::i,j,k
    real*8::heat_ZCIE,n(:),inTgas,Tgas
    real*8::v1,v2,v3,prev1,prev2,cH
    real*8::vv1_h,vv2_h,vv3_h,vv4_h,vv12_h,vv34_h,xGd
    real*8::x1(imax),x2(jmax),x3(kmax)
    real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
    real*8::v1min,v1max,v2min,v2max,v3min,v3max
    real*8,parameter::eps=1d-5

    Tgas = inTgas
    heat_ZCIE = 0d0

    !local copy of limits
    v1min = coolZCIEx1min
    v1max = coolZCIEx1max
    v2min = coolZCIEx2min
    v2max = coolZCIEx2max
    v3min = coolZCIEx3min
    v3max = coolZCIEx3max

    !local copy of variables arrays
    x1(:) = coolZCIEx1(:)
    x2(:) = coolZCIEx2(:)
    x3(:) = coolZCIEx3(:)

    ixd1(:) = coolZCIEixd1(:)
    ixd2(:) = coolZCIEixd2(:)
    ixd3(:) = coolZCIEixd3(:)

    !local variables
    cH = get_Hnuclei(n(:))

    !check if the abundance is close to zero to
    !avoid weird log evaluation
    if(cH.lt.1d-20)return

    v1 = Tgas           !Tgas
    v2 = cH             !total H number density
    v3 = phys_zredshift !redshift is linear

    !logs of variables
    v1 = log10(v1)
    v2 = log10(v2)

    !check limits
    if(v1>=v1max) v1 = v1max*(1d0-eps)
    if(v2>=v2max) v2 = v2max*(1d0-eps)
    if(v3>=v3max) v3 = v3max*(1d0-eps)

    if(v1<v1min) return
    if(v2<v2min) return
    if(v3<v3min) return

    !gets position of variable in the array
    i = (v1-v1min)*coolZCIEdvn1+1
    j = (v2-v2min)*coolZCIEdvn2+1
    k = (v3-v3min)*coolZCIEdvn3+1

    !precompute shared variables
    prev1 = (v1-x1(i))*ixd1(i)
    prev2 = (v2-x2(j))*ixd2(j)

    !linear interpolation on x1 for x2,x3
    vv1_h = prev1 * (heatZCIEy(k,j,i+1) - &
        heatZCIEy(k,j,i)) + heatZCIEy(k,j,i)
    !linear interpolation on x1 for x2+dx2,x3
    vv2_h = prev1 * (heatZCIEy(k,j+1,i+1) - &
        heatZCIEy(k,j+1,i)) + heatZCIEy(k,j+1,i)
    !linear interpolation on x2 for x3
    vv12_h = prev2 * (vv2_h - vv1_h) + vv1_h

    !linear interpolation on x1 for x2,x3+dx3
    vv3_h = prev1 * (heatZCIEy(k+1,j,i+1) - &
        heatZCIEy(k+1,j,i)) + heatZCIEy(k+1,j,i)
    !linear interpolation on x1 for x2+dx2,x3+dx3
    vv4_h = prev1 * (heatZCIEy(k+1,j+1,i+1) - &
        heatZCIEy(k+1,j+1,i)) + heatZCIEy(k+1,j+1,i)
    !linear interpolation on x2 for x3+dx3
    vv34_h = prev2 * (vv4_h - vv3_h) + vv3_h

    !linear interpolation on x3
    xGd = (v3-x3(k))*ixd3(k)*(vv34_h - &
        vv12_h) + vv12_h

    !Z cooling in erg/s/cm3
    heat_ZCIE = 1d1**xGd * cH * cH * total_Z

  end function heat_ZCIE
#ENDIFKROME


#IFKROME_useHeatingVisc
  !*************************
  !heating from viscosity (erg/s/cm3)
  ! requires user_nu (kinematic viscosity) and
  ! user_omega (keplerian orbital frequency)
  ! both from krome_user_commons
  function heat_visc(n,Tgas)
    use krome_commons
    use krome_user_commons
    use krome_subs
    implicit none
    real*8::n(:),Tgas,heat_visc
    real*8::m(nspec),rhogas

    n(idx_Tgas) = Tgas
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)

    heat_visc = 9d0/4d0 * user_nu * rhogas * user_omega * user_omega

  end function heat_visc
#ENDIFKROME


#IFKROME_useHeatingXRay
  !*************************
  !heating from xrays in erg/s/cm3
  function heat_XRay(n,Tgas,k)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_getphys
    use krome_fit
    implicit none
    real*8::n(:),Tgas,heat_Xray,k(:),ntot
    real*8::xheat_H,xheat_He,logH,logHe
    real*8::xe,ratexH,ratexHe,ncolH,ncolHe

    ntot = get_Hnuclei(n(:))
    xe = min(n(idx_e)/ntot,1d0)
    n(idx_Tgas) = Tgas
    !prepares logs for xrays
    ncolH = num2col(n(idx_H),n(:))
    ncolHe = num2col(n(idx_He),n(:))
    logH = log10(ncolH)
    logHe = log10(ncolHe)

    heat_Xray = 0d0
    xheat_H = fit_anytab2D(user_xheat_H_anytabx(:), &
         user_xheat_H_anytaby(:), &
         user_xheat_H_anytabz(:,:), &
         user_xheat_H_anytabxmul, &
         user_xheat_H_anytabymul, &
         logH,logHe-logH)
    xheat_He = fit_anytab2D(user_xheat_He_anytabx(:), &
         user_xheat_He_anytaby(:), &
         user_xheat_He_anytabz(:,:), &
         user_xheat_He_anytabxmul, &
         user_xheat_He_anytabymul, &
         logH,logHe-logH)

    !prepares varibles for xray photochemistry
    ratexH = 1d1**xheat_H * J21xray
    ratexHe = 1d1**xheat_He * J21xray

    heat_Xray = ratexH * n(idx_H)
    heat_Xray = heat_Xray + ratexHe * n(idx_He)
    heat_Xray = heat_Xray * .9971d0 * (1d0-(1d0-xe**.2663)**1.3163)

    heat_Xray = eV_to_erg * heat_Xray

  end function heat_XRay
#ENDIFKROME

#IFKROME_useHeatingPhotoDust
  !***************************
  function heat_photoDust(n,Tgas)
    !photoelectric effect from dust in erg/s/cm3
    !see Bakes&Tielens 1994 with a slight modification of Wolfire 2003
    !on the amount of absorbed ultraviolet energy.
    !This is for the local interstellar Habing flux and
    !without considering the recombination (which at this
    !radiation flux is indeed negligible)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8::heat_photoDust,n(:),Tgas,ntot,eps
    real*8::Ghab,z,psi

    ntot = get_Hnuclei(n(:))
    Ghab = 1.69d0 !habing flux, 1.69 is Draine78
#KROME_GhabG0
#KROME_GhabAv
    if(n(idx_e)>0d0) then
       psi = Ghab * sqrt(Tgas) / n(idx_e)
    else
       psi = 1d99
    end if
    eps = 4.9d-2 / (1d0 + 4d-3 * psi**.73) + &
         3.7d-2 * (Tgas * 1d-4)**.7 / (1d0 + 2d-4 * psi)
    z = #KROME_photoDustZ !metallicty wrt solar
    heat_photoDust = 1.3d-24*eps*Ghab*ntot*z

  end function heat_photoDust
#ENDIFKROME

#IFKROME_useHeatingPhotoDustNet
  !***************************
  function heat_netPhotoDust(n,Tgas)
    !photoelectric effect from dust in erg/s/cm3
    !including the recombination cooling and a generic radiation flux
    !eq. 42 and 44 in Bakes&Tielens, 1994
    ! dust2gas_ratio is D/D_sol, default assumes D/D_sol = Z/Z_sol
    use krome_commons
    use krome_subs
    use krome_constants
    use krome_getphys
    implicit none
    integer::i
    real*8::heat_netPhotoDust,n(:),Tgas,ntot,eps
    real*8::psi,recomb_cool,bet

    ntot = get_Hnuclei(n(:))
    bet = 0.735d0*(Tgas)**(-0.068)

    if(n(idx_e)>0d0) then
       psi = GHabing * sqrt(Tgas) / n(idx_e)
    else
       psi = 0d0
    end if

    !grains recombination cooling
    recomb_cool = 4.65d-30*Tgas**0.94*psi**bet &
         * n(idx_e)*n(idx_H)

    eps = 4.9d-2 / (1d0 + 4d-3 * psi**.73) + &
         3.7d-2 * (Tgas * 1d-4)**.7 / (1d0 + 2d-4 * psi)

    !net photoelectric heating
    heat_netPhotoDust = (1.3d-24*eps*GHabing*ntot-recomb_cool)*dust2gas_ratio

  end function heat_netPhotoDust
#ENDIFKROME

#IFKROME_useHeatingPhotoAv
  !******************************
  function heat_photoAv(n,Tgas,k)
    !heating from photoreactions using rate approximation (erg/s/cm3)
    use krome_commons
    use krome_user_commons
    use krome_subs
    use krome_getphys
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
    real*8::logH2,QH2,QH,QHe,ev2erg

    ev2erg = 1.60217662d-12
    Hfact = 2d1*ev2erg !erg

    !precompute log10(H2)
    logH2 = log10(max(n(idx_H2),1d-40))

    !init heating
    heat_CR = 0d0

    !heating per H ionization (eV)
    QH = 4.3d0 * ev2erg

    !heating per He ionization, same as H following Glassgold+2012
    QHe = QH

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
  function photo_heating(n,f2)
    !photo heating in erg/cm3/s using bin-based
    ! approach. Terms are computed in the
    ! krome_photo module
    use krome_commons
    use krome_constants
    implicit none
    real*8::photo_heating,n(:),f2

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
    use krome_getphys
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

    dd = sum(n(1:nmols)) !total number density

    !COMPRESSIONAL HEATING
    heat_compress = dd * boltzmann_erg * Tgas / user_tff !erg/s/cm3

  end function heat_compress
#ENDIFKROME

end module KROME_heating
