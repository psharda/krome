module krome_main
contains

#KROME_header

  !********************************
  !KROME main (interface to the solver library)
#IFKROME_useX
  subroutine krome(x,rhogas,Tgas,dt#KROME_dust_arguments)
#ELSEKROME
  subroutine krome(x,Tgas,dt#KROME_dust_arguments)
#ENDIFKROME
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    real*8::dt,x(nmols),rhogas,Tgas,mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,xdust(ndust)
    integer::icount,i
    
    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
#KROME_iwork_array
    real*8::atol(nspec),rtol(1)
#KROME_rwork_array
    logical::got_error


    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0.d0
    itol = 2 !both tolerances are scalar
    rtol = 1d-2 !relative tolerance (default: 1d-4)
    atol(:) = 1d-6 !absolute tolerance (default: 1d-40)
    itask = 1
    iopt = 0
    !MF=
    !  = 222 internal-generated JAC and sparsity
    !  = 121 user-provided JAC and internal generated sparsity
    !  =  22 internal-generated JAC but sparsity user-provided
    !  =  21 user-provided JAC and sparsity
#KROME_MF
    !end init DLSODES
    !****************************

#KROME_init_JAC
#KROME_init_IAC
    
    n(:) = 0.d0 !initialize densities
    
#IFKROME_useX
    mass(:) = get_mass() !get masses
    xin = sum(x) !store initial fractions
    !compute densities from fractions
    do i = 1,nspec
       if(mass(i)>0.d0) n(i) = rhogas * x(i) / mass(i)
    end do
#ELSEKROME
    n(1:nmols) = x(:)
#ENDIFKROME

    n(idx_Tgas) = Tgas !put temperature in the input array
    nold(:) = n(:) !store initial densities (needed for Jacobian)
    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

#IFKROME_useFlux
    rrmax = fex_check(n(:),Tgas)
    call flux_reduction(rrmax)
#ENDIFKROME

#IFKROME_check_mass_conservation
    mass(:) = get_mass() !get masses
    totmass = sum(n(:) * mass(:)) !calculate total mass
#ENDIFKROME 

    do
       icount = icount + 1
       !solve ODE
       CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
#IFKROME_report
       call krome_dump(n(:), rwork(:), iwork(:))
#ENDIFKROME
       !check DLSODES exit status
       if(istate==2) then
          exit !sucsessful integration
       elseif(istate==-1) then
          istate = 1 !exceeded internal max iterations
       elseif(istate==-5) then
          istate = 3 !wrong sparsity
       else
          call XSETF(1)!turn on verbosity
          if(got_error) then
             print *,"ERROR: wrong solver exit status!"
             print *,"istate:",istate
             print *,"SEE KROME_ERROR_REPORT file"
             call krome_dump(n(:), rwork(:), iwork(:))
             stop
          end if
          got_error = .true.
          istate = 1
       end if
    end do

#IFKROME_check_mass_conservation
    if(abs(1.d0-totmass/sum(n(:) * mass(:)))>1d-3) then
       print *,"ERROR: mass conservation failure!"
       print *, tloc,totmass,sum(n(:) * mass(:))
    end if
#ENDIFKROME

    !avoid negative species
    do i=1,nspec
       n(i) = max(n(i),0.d0)
    end do

#IFKROME_useX
    x(:) = mass(1:nmols)*n(1:nmols)/rhogas !return to fractions
    x(:) = x(:) / sum(x) * xin !force mass conservation
#ELSEKROME
    x(:) = n(1:nmols)
#ENDIFKROME

    Tgas = n(idx_Tgas) !get new temperature
  end subroutine krome

  !*******************************
  subroutine krome_dump(n,rwork,iwork)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_reduction
    integer::fnum,i,iwork(:)
    real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc
    character*16::names(nspec),FMTi,FMTr
    fnum = 99
    open(fnum,FILE="KROME_ERROR_REPORT",status="replace")
    
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
    k(:) = coe_tab(n(:))
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
    rrmax = fex_check(n(:), n(idx_Tgas))
    write(fnum,*)
    write(fnum,*) "Reaction magnitude [k*n1*n2*n3]"
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12)') "#","flux","flux %"
    write(fnum,*) "**********************"    
    do i=1,nrea
       rperc = 0.d0
       if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
       write(fnum,'(I5,2E12.3e3)') i,arr_flux(i),rperc
    end do
    write(fnum,*) "**********************"
    write(fnum,*)

    !SOLVER
    FMTr = "(a30,E16.7e3)"
    FMTi = "(a30,I10)"
    write(fnum,*) "Solver-related information:"
    write(fnum,FMTr) "step size last",rwork(11)
    write(fnum,FMTr) "step size attempt",rwork(12)
    write(fnum,FMTr) "time current",rwork(13)
    write(fnum,FMTr) "tol scale factor",rwork(14)
    write(fnum,FMTi) "numeber of steps",iwork(11)
    write(fnum,FMTi) "call to fex",iwork(12)
    write(fnum,FMTi) "call to jex",iwork(13)
    write(fnum,FMTi) "last order used",iwork(14)
    write(fnum,FMTi) "order attempt",iwork(15)
    write(fnum,FMTi) "idx largest error",iwork(16)
    write(fnum,FMTi) "RWORK size required",iwork(17)
    write(fnum,FMTi) "IWORK size required",iwork(18)
    write(fnum,FMTi) "NNZ in Jac",iwork(19)
    write(fnum,FMTi) "extra fex to compute jac",iwork(20)
    write(fnum,FMTi) "number of LU decomp",iwork(21)
    write(fnum,FMTi) "base address in RWORK",iwork(22)
    write(fnum,FMTi) "base address of IAN",iwork(23)
    write(fnum,FMTi) "base address of JAN",iwork(24)
    write(fnum,FMTi) "NNZ in lower LU",iwork(25)
    write(fnum,FMTi) "NNZ in upper LU",iwork(21)    
    write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
    write(fnum,*) 
    write(fnum,*) "END KROME ERROR REPORT"
    write(fnum,*)
    close(fnum)
  end subroutine krome_dump

  !********************************
  subroutine krome_init()
    use krome_tabs
    use krome_subs
    use krome_reduction

    call load_arrays

#IFKROME_useTabs
    call make_ktab()
    call check_tabs()
#ENDIFKROME
#IFKROME_useTopology
    call krome_authub()
#ENDIFKROME

  end subroutine krome_init

  !****************************
  function krome_get_coe(n)
    !krome_get_coe: public interface to obtain rate coefficients
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
    real*8::krome_get_coe(nrea),n(:)
    krome_get_coe(:) = coe_tab(n(:))
  end function krome_get_coe

  !****************************
  function krome_get_coeT(Tgas)
    !krome_get_coeT: public interface to obtain rate coefficients
    ! with argument Tgas only
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
    real*8::krome_get_coeT(nrea),n(nspec),Tgas
    n(idx_Tgas) = Tgas
    krome_get_coeT(:) = coe_tab(n(:))
  end function krome_get_coeT

end module krome_main
