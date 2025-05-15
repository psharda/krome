module krome_ode
contains

#KROME_header

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
#IFKROME_useDust
    use krome_dust
#ENDIFKROME
#IFKROME_useCoolingDustSemenov
    use krome_dust
#ENDIFKROME
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq),nH
#KROME_iceODEVariables
#KROME_dustSumVariables
#KROME_implicit_variables
#KROME_flux_variables
#KROME_initcoevars

#KROME_coevars

    n(:) = nin(:)
    ntot = sum(n(1:nmols))
#IFKROME_popsicle_ice
    ntot = sum(n(1:nmols)) - n(idx_CO_total) - n(idx_H2O_total)
#ENDIFKROME_popsicle_ice
    nH = get_Hnuclei(n(:))
    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

#KROME_compute_electrons

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),phys_Tcmb)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

#IFKROME_tigressNCR
    !If Eq C chemistry used, update final abundances for these to the equilibrium one
    call steadystate_tigressNCR(n(:), Tgas)
#ENDIFKROME_tigressNCR

#IFKROME_shieldHabingDust
    call calcHabingThick(n(:),Tgas)
#ENDIFKROME

#KROME_Tdust_limits

    k(:) = coe_tab(n(:)) !compute coefficients

#KROME_iceODEDefinitions

#KROME_H2pdRate

#KROME_photobins_compute_thick

#IFKROME_useDust
    vgas = sqrt(kvgas_erg*Tgas)
    ntot = sum(n(1:nmols))

    #KROME_calc_Tdust

#ENDIFKROME

#IFKROME_useCoolingDustSemenov
    call compute_Semenov_Tdust(n(:),Tgas)
#ENDIFKROME

#KROME_dust_H2

#KROME_ODE

#IFKROME_use_thermo_toggle
    if(krome_thermo_toggle>0) then
#ENDIFKROME


#IFKROME_applyElementConservation_popsicle_semenov
    !No Photo
    !species: E, H-, D-, H, He, H2, HD, D, O, OH, H2O, O2, CH, CO, C, CO2, CH2, CH3, CH4, H+, He+, H2+, D+, HD+, O+, OH+, H2O+, H3O+, O2+, C+, CO+, He++

    !H conservation
    dn(idx_H) = -1d0 * (dn(idx_Hk) + 2*dn(idx_H2) + dn(idx_HD) + dn(idx_OH) + 2*dn(idx_H2O_total) + dn(idx_CH) + 2*dn(idx_CH2) + 3*dn(idx_CH3) + 4*dn(idx_CH4) + &
                        dn(idx_Hj) + 2*dn(idx_H2j) + dn(idx_HDj) + dn(idx_OHj) + 2*dn(idx_H2Oj) + 3*dn(idx_H3Oj))

    !He conservation
    dn(idx_He) = -1d0 * (dn(idx_HEj) + dn(idx_Hejj))

    !C conservation
    dn(idx_C) = -1d0 * (dn(idx_CH) + dn(idx_CO_total) + dn(idx_CO2) + dn(idx_CH2) + dn(idx_CH3) + dn(idx_CH4) + dn(idx_Cj) + dn(idx_COj))

    !O conservation
    dn(idx_O) = -1d0 * (dn(idx_OH) + dn(idx_H2O_total) + 2*dn(idx_O2) + dn(idx_CO_total) + 2*dn(idx_CO2) + dn(idx_Oj) + dn(idx_OHj) + dn(idx_H2Oj) + dn(idx_H3Oj) + &
                        2*dn(idx_O2j) + dn(idx_COj))

    !D conservation
    dn(idx_D) = -1d0 * (dn(idx_Dk) + dn(idx_Dj) + dn(idx_HD) + dn(idx_HDj))
#ENDIFKROME

#IFKROME_applyElementConservation_popsicle_semenov_photo_full
    !With Photo (Full Network)
    !species: E, H-, D-, H, He, H2, HD, D, O, OH, H2O, O2, CH, CO, C, CO2, CH2, CH3, CH4, H+, He+, H2+, D+, HD+, O+, OH+, H2O+, H3O+, O2+, C+, CO+, He++, CH+, CH2+, H3+

    !H conservation
    dn(idx_H) = -1d0 * (dn(idx_Hk) + 2*dn(idx_H2) + dn(idx_HD) + dn(idx_OH) + 2*dn(idx_H2O_total) + dn(idx_CH) + 2*dn(idx_CH2) + 3*dn(idx_CH3) + 4*dn(idx_CH4) + &
                        dn(idx_Hj) + 2*dn(idx_H2j) + dn(idx_HDj) + dn(idx_OHj) + 2*dn(idx_H2Oj) + 3*dn(idx_H3Oj) + 3*dn(idx_H3j) + dn(idx_CHj) + 2*dn(idx_CH2j))

    !He conservation
    dn(idx_He) = -1d0 * (dn(idx_HEj) + dn(idx_Hejj))

    !C conservation
    dn(idx_C) = -1d0 * (dn(idx_CH) + dn(idx_CO_total) + dn(idx_CO2) + dn(idx_CH2) + dn(idx_CH3) + dn(idx_CH4) + dn(idx_Cj) + dn(idx_COj) + dn(idx_CHj) + dn(idx_CH2j))

    !O conservation
    dn(idx_O) = -1d0 * (dn(idx_OH) + dn(idx_H2O_total) + 2*dn(idx_O2) + dn(idx_CO_total) + 2*dn(idx_CO2) + dn(idx_Oj) + dn(idx_OHj) + dn(idx_H2Oj) + dn(idx_H3Oj) + &
                        2*dn(idx_O2j) + dn(idx_COj))

    !D conservation
    dn(idx_D) = -1d0 * (dn(idx_Dk) + dn(idx_Dj) + dn(idx_HD) + dn(idx_HDj))
#ENDIFKROME


#IFKROME_applyElementConservation_popsicle_semenov_photo_gow
    !species: E, H, H2, OH, CH, H+, H2+, H3+, HCO+, He, He+, He++, C, CH, CO, C+, O, O+, Si, Si+

    !H conservation
    dn(idx_H) = -1d0 * (2*dn(idx_H2) + dn(idx_OH) + dn(idx_CH) + &
                        dn(idx_Hj) + 2*dn(idx_H2j) + 3*dn(idx_H3j) + dn(idx_HCOj))

    !He conservation
    dn(idx_He) = -1d0 * (dn(idx_HEj) + dn(idx_Hejj))

    !C conservation
    dn(idx_C) = -1d0 * (dn(idx_CH) + dn(idx_CO) + dn(idx_Cj) + dn(idx_HCOj))

    !O conservation
    dn(idx_O) = -1d0 * (dn(idx_OH) + dn(idx_CO) + dn(idx_Oj) + &
                        dn(idx_HCOj))

    !Si conservation
    dn(idx_SI) = -1d0 * dn(idx_SIj)

#ENDIFKROME


#IFKROME_use_thermo
       krome_gamma = gamma_index(n(:))

       dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
            - cooling(n(:), Tgas) #KROME_coolingQuench #KROME_coolfloor) &
            * (krome_gamma - 1.d0) / boltzmann_erg / ntot
#ENDIFKROME

#IFKROME_use_thermo_toggle
    end if
#ENDIFKROME


#IFKROME_usedTdust
    dn(nmols+ndust+1:nmols+2*ndust) = get_dTdust(n(:),dn(idx_Tgas),vgas,ntot)
#ENDIFKROME

#KROME_ODEModifier

#KROME_odeConstant

#IFKROME_report
       call krome_ode_dump(n(:), k(:))
       write(98,'(999E12.3e3)') tt,n(:)
       write(97,'(999E12.3e3)') tt,dn(:)
#ENDIFKROME

       last_coe(:) = k(:)

  end subroutine fex

#IFKROME_tigressNCR
  !Return the steady state abundances of C, C+, O, O+, CO given CR/photo rates and non-eq H abundances
  !Refer to Kim et al. 2023, ApJS, 264, 1 for more details
  subroutine steadystate_tigressNCR(n,Tgas)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    real*8::n(:), Tgas
    real*8::ne,xCtot,xOtot, pion_C, crion_C, alpha17, beta17, gamma17, krr, kdr1, kdr2, kdr3, kdr4, kdr5, &
            kdr, nCplus, kCplusH2, nOplus, nCO, nCO_crit, nH, ntot
    real, parameter :: dissCO_draine = 2.592d-10 !Base CO dissociation rate for Draine ISRF field

    ntot = sum(n(1:nmols))
    nH = get_Hnuclei(n(:))

    !Set steady state C/C+/O/O+ abundances
    xCtot = 1.6d-4 * phys_metallicity !Total C nuclei abundance w.r.t H nuclei (NOTE: Hardcoded for now)
    xOtot = 3.2d-4 * phys_metallicity !Total O nuclei abundance w.r.t H nuclei (NOTE: Hardcoded for now)
    !Definitions for the radiative+dielectronic recombinations of C+
    alpha17=sqrt(Tgas/6.67d-3)
    beta17=sqrt(Tgas/1.943d6)
    gamma17=0.7849d0 + 0.1597*exp(-49550d0/Tgas)
    !Radiative recombination rate of C+
    krr = 2.995d-9 / (alpha17*((1d0+alpha17)**(1d0-gamma17))*(1d0+beta17)**(1d0+gamma17))
    kdr1 = 6.346d-9*exp(-12.17d0/Tgas)
    kdr2 = 9.793d-9*exp(-73.80d0/Tgas)
    kdr3 = 1.634d-6*exp(-15230d0/Tgas)
    kdr4 = 8.369d-4*exp(-1.207d5/Tgas)
    kdr5 = 3.355d-4*exp(-2.144d5/Tgas)
    !Dielectronic recombination rate of C+
    kdr = (Tgas**(-1.5d0))*(kdr1+kdr2+kdr3+kdr4+kdr5)
    ne = n(idx_E)
    !Photoionization rate of C (Eq. 23; Kim+ 2023)
    pion_C = user_ionC + 520*2*(n(idx_H2)/nH)*user_crate !Photoionization + CR ionization from Heays et al. 2017
    !Cosmic-ray ionization rate of C
    crion_C = 3.85 * user_crate !cosmic-ray ionization rate of C
    !Reaction rate of the charge transfer reaction C+ + H2-> CH2+ (reaction 9 in GOW)
    kCplusH2 = 2.31d-13*(Tgas)**(-1.3)*exp(-23/Tgas)

    !Equation 24 in Kim+2023
    nCplus = xCtot * nH * (pion_C + crion_C)/(pion_C + crion_C+ntot*C_recombination_on_dust(n,Tgas) &
        + (krr+kdr)*ne + kCplusH2*n(idx_H2))

    !O/O+ steady state solution is simply equal to the H/H+ ratio due to efficient charge transfer
    nOplus = xOtot * nH * n(idx_Hj)/nH

    !Critical density above which xCO/xC,tot > 0.5 ; Equation 26 in Kim+2023
    nCO_crit = (4d3 * dust2gas_ratio * (user_crate/1d-16)**(-2d0)) ** ((user_chiCO)**(1d0/3d0)) * &
                (5d1 * (user_crate/1d-16))/(dust2gas_ratio**(1.4d0))
    !Equation 25 in Kim+2023
    nCO = xCtot * nH * 2d0 * (n(idx_H2)/nH) * (1d0 - max(nCplus/(xCtot*nH),nOplus/(xOtot*nH)))/(1d0 + (nCO_crit/ntot)**2d0)

    n(idx_Cj) = max(nCplus,0d0)
    n(idx_CO) = max(nCO,0d0)
    !Obtain nC from closure (only C, C+ and CO are considered)
    n(idx_C) = max(xCtot * nH - nCplus - nCO,0d0)

    nOplus = max(nOplus,0d0)
    !Obtain nO from closure (only O and O+ are considered)
    n(idx_O) = max(xOtot * nH - nOplus -nCO,0d0)
    
    return
  end subroutine steadystate_tigressNCR
#ENDIFKROME_tigressNCR

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

#KROME_compute_electrons

#IFKROME_use_thermo
    krome_gamma = gamma_index(n(:))
#ENDIFKROME

    k(:) = last_coe(:) !get rate coefficients

#KROME_JAC_PD

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

#KROME_JAC_PDX

  end subroutine jex

#IFKROME_report

  !*******************************
  subroutine krome_ode_dump(n,k)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_getphys
    integer::fnum,i
    real*8::n(:),rrmax,k(:),kmax,kperc
    character*16::names(nspec)
    character*50::rnames(nrea)
    fnum = 99
    open(fnum,FILE="KROME_ODE_REPORT",status="replace")

    names(:) = get_names()

    write(fnum,*) "KROME ERROR REPORT"
    write(fnum,*)
    !SPECIES
    write(fnum,*) "Species aboundances"
    write(fnum,*) "**********************"
    write(fnum,'(a5,a20,a12)') "#","name","qty"
    write(fnum,*) "**********************"
    do i=1,nspec
       write(fnum,'(I5,a20,E12.3e3)') i,names(i),n(i)
    end do
    write(fnum,*) "**********************"

    !RATE COEFFIECIENTS
    kmax = maxval(k)
    write(fnum,*)
    write(fnum,*) "Rate coefficients at Tgas",n(idx_Tgas)
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12)') "#","k","k %"
    write(fnum,*) "**********************"
    do i=1,nrea
       kperc = 0.d0
       if(kmax>0.d0) kperc = k(i)*1d2/kmax
       write(fnum,'(I5,2E12.3e3)') i,k(i),kperc
    end do
    write(fnum,*) "**********************"
    write(fnum,*)

    !FLUXES
    call load_arrays
    rnames(:) = get_rnames()
    write(fnum,*)
    write(fnum,*) "Reaction magnitude [k*n1*n2*n3]"
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12)') "#","flux"
    write(fnum,*) "**********************"
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
#KROME_report_flux
    end do
    write(fnum,*) "**********************"
    write(fnum,*)

    write(fnum,*) "END KROME ODE REPORT"
    write(fnum,*)
    close(fnum)
  end subroutine krome_ode_dump
#ENDIFKROME

end module krome_ode
