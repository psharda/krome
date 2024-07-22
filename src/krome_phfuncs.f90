!This module mainly contains shielding routine and
! function to initialize radiation background (e.g. Planck).
module krome_phfuncs
contains

#KROME_header

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
#IFKROME_hasHI
    NHtot  = NHtot + num2col(n(idx_H),n(:))
#ENDIFKROME
#IFKROME_hasHII
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
#ENDIFKROME
#IFKROME_hasH2I
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))
#ENDIFKROME

    shield_dust = exp(-sigma_d*NHtot)

  end function shield_dust

#IFKROME_usePhotoBins

  !*******************
  !apply a shielding to Habing flux
  subroutine calcHabingThick(n,Tgas)
    use krome_commons
    implicit none
    real*8::getHabingThick,n(:),Tgas

    GHabing = GHabing_thin * shield_dust(n(:),Tgas,0.665d0)

  end subroutine calcHabingThick

  !*********************
  !return the ratio between the current flux an Draine's
  function get_ratioFluxDraine()
    implicit none
    real*8::get_ratioFluxDraine

    !7.95d-8 eV/cm2/sr is the integrated Draine flux
    get_ratioFluxDraine = get_integratedFlux()/7.95d-8

  end function get_ratioFluxDraine

  !**********************
  !return the curred integrated flux (eV/cm2/sr)
  ! as I(E)/E*dE
  function get_integratedFlux()
    use krome_commons
    implicit none
    integer::j
    real*8::get_integratedFlux,dE

    get_integratedFlux = 0d0
    do j=1,nPhotoBins
       dE = photoBinEdelta(j)
       get_integratedFlux = get_integratedFlux &
            + photoBinJ(j)*dE/photoBinEmid(j)
    end do

  end function get_integratedFlux

#ENDIFKROME

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
       planckBB = 2d0*x**3/planck_eV**2/clight**2 &
            / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
       a = 2d0/planck_eV**2/clight**2
       planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

#IFKROME_useShieldingDB96
    !compute shielding from Draine+Bertoldi 1996
    krome_fshield = calc_H2shieldDB96(n(:), Tgas)
#ENDIFKROME

#IFKROME_useShieldingWG11
    !compute shielding from Wolcott+Greene 2011
    krome_fshield =  calc_H2shieldWG11(n(:), Tgas)
#ENDIFKROME

#IFKROME_useShieldingR14
    !compute shielding from Richings+ 2014
    krome_fshield =  calc_H2shieldR14(n(:), Tgas)
#ENDIFKROME

  end function krome_fshield

  !***********************
  !shielding function for CO selected with -shielding_CO option
  function krome_fshield_CO(n,Tgas)
    implicit none
    real*8::krome_fshield_CO,n(:),Tgas

    krome_fshield_CO = 1d0 !default shielding value

#IFKROME_useShieldingCO
    !compute shielding from Visser, van Dishoeck and Black 2009
    krome_fshield_CO = calc_COshieldVvDB09(n(:), Tgas)
#ENDIFKROME

  end function krome_fshield_CO

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
       fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
       fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
         0.035d0/sqrt(1d0+x) * &
         exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

#IFKROME_useShieldingDB96
  !************************
  !calculate the self-shielding factor, following Draine&Bertoldi 1996
  !NOTE: this function is suited for collapse. Use with caution!
  function calc_H2shieldDB96(n,Tgas)
    use krome_commons
    real*8::n(nspec),Tgas,calc_H2shieldDB96,N_H2, nH2

    !check on H2 abundances to avoid
    ! weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))

!    N_H2 = nH2*get_jeans_length(n(:),Tgas)*0.5d0  !column density (cm-2)
    N_H2  =  2d0 * num2col(nH2,n(:))

    calc_H2shieldDB96 = min(1.d0, (N_H2*1.d-14)**(-0.75d0))

  end function calc_H2shieldDB96
#ENDIFKROME

#IFKROME_useShieldingWG11
  !************************
  !calculate the self-shielding factor, following Wolcott&Greene 2011
  !NOTE: this function is suited for collapse. Use with caution!
  function calc_H2shieldWG11(n,Tgas)
    use krome_commons
    use krome_constants
    use krome_getphys
    real*8::n(nspec),Tgas,calc_H2shieldWG11,N_H2,nH2
    real*8::xN_H2,b5,H_mass

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))

    N_H2  =  2d0 * num2col(nH2,n(:))

!    N_H2 = nH2*get_jeans_length(n(:) ,Tgas)*0.5d0  !column density (cm-2)
    xN_H2 = N_H2*2d-15 !normalized column density (#), 2d-15=1/5d14
    H_mass = p_mass+e_mass !H mass in g

    !doppler broadening parameter b divided by 1d5 cm/s (#)
    b5 = ((boltzmann_erg*Tgas/H_mass)**0.5d0)*1.d-5
    calc_H2shieldWG11 = 0.965d0/(1.d0+xN_H2/b5)**1.1d0 &
         + (0.035d0/(1.d0+xN_H2)**0.5d0) &
         * exp(-8.5d-4*(1.d0+xN_H2)**0.5d0)

  end function calc_H2shieldWG11
#ENDIFKROME

#IFKROME_useShieldingR14
  !Temperature-dependent self-shielding as reported in Richings+2014.
  function calc_H2shieldR14(n,Tgas)
    use krome_commons
    use krome_constants
    use krome_getphys
    real*8::n(nspec),Tgas,calc_H2shieldR14,N_H2,nH2
    real*8::xN_H2,b5,H_mass,bturb,btherm2
    real*8::alpha,omegaH2,Ncrit

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))

    N_H2  =  2d0 * num2col(nH2,n(:))

!    N_H2 = nH2*get_jeans_length(n(:) ,Tgas)*0.5d0  !column density (cm-2)
    H_mass = p_mass+e_mass !H mass in g
    bturb = 7.1d0*km_to_cm !turbulent Doppler broadening parameter in cm/s
    btherm2 = boltzmann_erg*Tgas/H_mass !thermal Doppler broadening parameter cm/s

    !doppler broadening parameter b divided by 1d5 cm/s (#)
    b5 = ((btherm2 + bturb**2d0)**0.5)*1.d-5
    omegaH2 = 0.013d0*(1d0+(Tgas/2.7d3)**1.3)**(1.0/1.3)*exp(-(Tgas/3.9d3)**14.6)

    if(Tgas<3d3)then
      alpha = 1.4
      Ncrit = 1.3d0*(1d0+(Tgas/6d2)**0.8)
    elseif(Tgas>=3d3.or.Tgas<4d3)then
      alpha = (Tgas/4.5d3)**(-0.8)
      Ncrit = (Tgas/4.76d3)**(-3.8)
    else
      alpha = 1.1
      Ncrit = 2.d0
    endif

    xN_H2 = N_H2*1d-14/Ncrit !normalized column density (#)

    calc_H2shieldR14 = (1d0-omegaH2) / (1d0+xN_H2/b5)**alpha &
         * exp(-5d-7*(1d0+xN_H2)) &
         + (omegaH2/sqrt(1d0+xN_H2)) * exp(-8.5d-4*sqrt(1d0+xN_H2))

  end function calc_H2shieldR14
#ENDIFKROME

#IFKROME_useShieldingCO
  !************************
  !calculate the shielding factor for CO, following Visser, van Dishoeck and Black 2009
  !Shielding by CO and H2. Does not include dust
  function calc_COshieldVvDB09(n,Tgas)
    use krome_commons
    use krome_constants
    use krome_getphys
    use krome_fit
    real*8::n(nspec),Tgas,calc_COshieldVvDB09,N_H2,nH2
    real*8::N_CO,ncO,clipped_x,clipped_y
    real*8 :: x(8), y(6), z(8, 6)

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))
    nCO = max(1d-40, n(idx_CO))

    N_H2  =  2d0 * num2col(nH2,n(:))
    N_CO = num2col(nCO,n(:))

    !Table 5, first section
    x = (/ 0d0, 13d0, 14d0, 15d0, 16d0, 17d0, 18d0, 19d0 /) !N_CO
    y = (/ 0d0, 19d0, 20d0, 21d0, 22d0, 23d0 /) !N_H2
    z(:, 1) = (/ 1d0, 8.080d-1, 5.250d-1, 2.434d-1, 5.467d-2, 1.362d-2, 3.378d-3, 5.240d-5 /)
    z(:, 2) = (/ 8.176d-1, 6.347d-1, 3.891d-1, 1.787d-1, 4.297d-2, 1.152d-2, 2.922d-3, 4.662d-4 /)
    z(:, 3) = (/ 7.223d-1, 5.624d-1, 3.434d-1, 1.540d-1, 3.515d-2, 9.231d-3, 2.388d-3, 3.899d-4 /)
    z(:, 4) = (/ 3.260d-1, 2.810d-1, 1.953d-1, 8.726d-2, 1.907d-2, 4.768d-3, 1.150d-3, 1.941d-4 /)
    z(:, 5) = (/ 1.108d-2, 1.081d-2, 9.033d-3, 4.441d-3, 1.102d-3, 2.644d-4, 7.329d-5, 1.437d-5 /)
    z(:, 6) = (/ 3.938d-7, 3.938d-7, 3.936d-7, 3.923d-7, 3.901d-7, 3.893d-7, 3.890d-7, 3.875d-7 /)


    !Clip logNCO and logNH2 to the ranges in the data
    clipped_x = max(x(1), min(log10(N_CO), x(8)))
    clipped_y = max(y(1), min(log10(N_H2), y(6)))

    !Interpolate on log10 (z)
    calc_COshieldVvDB09 = 1d1**interpolate2D(x, y, log10(z), clipped_x, clipped_y)

  end function calc_COshieldVvDB09
#ENDIFKROME

end module krome_phfuncs
