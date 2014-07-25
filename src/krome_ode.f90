module krome_ode
contains

#KROME_header

  subroutine fex(neq,tt,n,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
#IFKROME_useDust
    use krome_dust
#ENDIFKROME
    implicit none
    integer::neq
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd
#KROME_dustSumVariables
#KROME_implicit_variables
#KROME_flux_variables
#KROME_initcoevars

#KROME_coevars

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    Tgas = max(n(idx_Tgas), 2.73d0) !get temperature
    
    k(:) = coe_tab(n(:)) !compute coefficients

#KROME_photobins_compute_thick

#IFKROME_useDust
    vgas = sqrt(kvgas_erg*Tgas)
    ntot = sum(n(1:nmols))

    #KROME_calc_Tdust

    #KROME_dust_H2
#ENDIFKROME

#KROME_ODE

#IFKROME_use_thermo_toggle
    if(krome_thermo_toggle>0) then
#ENDIFKROME

#IFKROME_use_thermo
       krome_gamma = gamma_index(n(:))
       
       dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
            - cooling(n(:), Tgas)#KROME_cool_cmb_floor) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
#ENDIFKROME
       
#IFKROME_use_thermo_toggle
    end if
#ENDIFKROME

#KROME_ODEModifier

#KROME_odeConstant

#IFKROME_report
       call krome_ode_dump(n(:), k(:))
       write(98,'(999E12.3e3)') tt,n(:)
       write(97,'(999E12.3e3)') tt,dn(:)
#ENDIFKROME
       
       jac_dnold(:) = jac_dn(:) !store previous dn for explicit jacobian
       jac_nold(:) = n(:) !store previous n for explicit jacobian
       jac_dn(:) = dn(:) !store current n

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

#IFKROME_use_thermo
    krome_gamma = gamma_index(n(:))
#ENDIFKROME
    
    if(j.ne.idx_Tgas) k(:) = coe_tab(n(:)) !compute coefficients
    
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
