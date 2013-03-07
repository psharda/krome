module krome_ode
contains

#KROME_header

  subroutine fex(neq,tt,n,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_tabs
    implicit none
    integer::neq
    real*8::tt,dn(neq),n(neq),k(nrea)
    real*8::gamma,Tgas
#KROME_implicit_variables

    dn(:) = 0.d0 !initialize differentials
    Tgas = max(n(idx_Tgas), 2.73d0) !get temperature
    k(:) = coe_tab(n(:)) !compute coefficients

#KROME_ODE
    
#IFKROME_use_cooling
    gamma = 5.d0/3.d0
    dn(idx_Tgas) = (heating(n(:), Tgas) - cooling(n(:), Tgas)) &
         * (gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
#ENDIFKROME

#KROME_odeConstant

#IFKROME_report
       call krome_ode_dump(n(:), k(:))
       write(98,'(999E12.3e3)') tt,n(:)
       write(97,'(999E12.3e3)') tt,dn(:)
#ENDIFKROME

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    
    Tgas = n(idx_Tgas)
    
#KROME_JAC_PD
    
    return
  end subroutine jes

  !****************************
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
       write(fnum,'(I5,E12.3e3,a2,a50)') i,k(i)*n(arr_r1(i))*n(arr_r2(i))*n(arr_r3(i)),"",rnames(i)
    end do
    write(fnum,*) "**********************"
    write(fnum,*)

    write(fnum,*) "END KROME ODE REPORT"
    write(fnum,*)
    close(fnum)
  end subroutine krome_ode_dump

end module krome_ode
