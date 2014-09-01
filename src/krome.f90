module krome_main

  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

#KROME_header

  !********************************
  !KROME main (interface to the solver library)
#IFKROME_useX
  subroutine krome(x,rhogas,Tgas,dt)!#KROME_dust_arguments)
#ELSEKROME
  subroutine krome(x,Tgas,dt)!#KROME_dust_arguments)
#ENDIFKROME
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    real*8::dt,x(nmols),rhogas,Tgas,mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec)
    integer::icount,i,ierr,icount_max
    
    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
#KROME_iwork_array
    real*8::atol(nspec),rtol(nspec)
#KROME_rwork_array
    logical::got_error,equil


    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0.d0
    itol = 4 !both tolerances are scalar
    rtol(:) = #KROME_RTOL !relative tolerance
    atol(:) = #KROME_ATOL !absolute tolerance
    icount_max = 100 !maximum number of iterations

#KROME_custom_RTOL

#KROME_custom_ATOL

    itask = 1
    iopt = 0

#KROME_maxord

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

    ierr = 0 !error flag, zero==OK!
    n(:) = 0.d0 !initialize densities
    
#IFKROME_useX
    mass(:) = get_mass() !get masses
    xin = sum(x) !store initial fractions
    !compute densities from fractions
    do i = 1,nmols
       if(mass(i)>0.d0) n(i) = rhogas * x(i) / mass(i)
    end do
#ELSEKROME
    n(1:nmols) = x(:)
#ENDIFKROME

    n(idx_Tgas) = Tgas !put temperature in the input array
    
#IFKROME_useDust
    n(nmols+1:nmols+ndust) = xdust(:) !get dust abundances
#ENDIFKROME


    jac_nold(:) = n(:) !store initial densities (finite difference for Jacobian)
    jac_dn(:) = 0.d0
    jac_dnold(:) = 0.d0
    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

#IFKROME_check_mass_conservation
    mass(:) = get_mass() !get masses
    totmass = sum(n(:) * mass(:)) !calculate total mass
#ENDIFKROME

    !store initial values
    ni(:) = n(:)

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
       icount = icount + 1
       !solve ODE
       CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)
#IFKROME_report
       call krome_dump(n(:), rwork(:), iwork(:), ni(:))
#ENDIFKROME
       krome_call_to_fex = krome_call_to_fex + IWORK(12)
       !check DLSODES exit status
       if(istate==2) then
          exit !sucsessful integration
       elseif(istate==-1) then
          istate = 1 !exceeded internal max iterations
       elseif(istate==-5 .or. istate==-4) then
          istate = 3 !wrong sparsity recompute
#IFKROME_noierr
       else
          call XSETF(1)!turn on verbosity
          got_error = .true.
          istate = 1
       end if
       
       if(got_error.or.icount>icount_max) then
          if (krome_mpi_rank>0) then
            print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
            print *,krome_mpi_rank,"istate:",istate
            print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
          else
            print *,"ERROR: wrong solver exit status!"
            print *,"istate:",istate
            print *,"SEE KROME_ERROR_REPORT file"
          end if
          call krome_dump(n(:), rwork(:), iwork(:), ni(:))
          stop
       end if
#ENDIFKROME
#IFKROME_ierr
       else
          !store the istate in ierr and exit from the loop
          ierr = istate
          exit
       end if
#ENDIFKROME

#IFKROME_useEquilibrium
       !try to determine if the system has reached a steady equilibrium
       equil = .true.
       do i=1,nspec
          if(n(i)>1d-40) then
             if(abs(n(i)-n_old(i))/n(i)>1d-6) then
                equil = .false.
                exit
             end if
          end if
       end do
       
       if(equil) exit
       n_old(:) = n(:)
#ENDIFKROME
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

#IFKROME_conserve
    n(:) = conserve(n(:),ni(:))
#ENDIFKROME

#IFKROME_useX
    x(:) = mass(1:nmols)*n(1:nmols)/rhogas !return to fractions
    x(:) = x(:) / sum(x) * xin !force mass conservation
#ELSEKROME
    !returns to user array
    x(:) = n(1:nmols)
#ENDIFKROME

#IFKROME_useDust
    xdust(:) = n(nmols+1:nmols+ndust)
#ENDIFKROME

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome


  !*********************************
  !integrates to equilibrium using constant temperature
  subroutine krome_equilibrium(x,Tgas)
    use krome_commons
    use krome_constants
    implicit none
    integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
    integer::i,imax
    real*8::tloc,x(:),Tgas,n(nspec),dt
#KROME_iwork_array
    real*8::atol(nspec),rtol(nspec)
#KROME_rwork_array

    call XSETF(0)!toggle solver verbosity
    meth = 2
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0.d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1d-6 !relative tolerance
    atol(:) = 1d-10 !absolute tolerance

    !for DLSODES options see its manual
    iopt = 0
    itask = 1
    istate = 1

    mf = 222 !internally evaluated sparsity and jacobian
    tloc = 0d0 !initial time
    dt = seconds_per_year * 1d8
    
    !copy into array
    n(nmols+1:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    
    imax = 1000

    do i=1,imax
       !solve ODE
       CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
       if(istate==2) then
          exit
       else
          istate=1
       end if
    end do
    !check errors
    if(istate.ne.2) then
       print *,"ERROR: no equilibrium found!"
       stop
    end if
    x(:) = n(1:nmols)

  end subroutine krome_equilibrium

  !********************
  !dummy jacobian
  subroutine jcn_dummy()
    implicit none
  end subroutine jcn_dummy

  !*******************
  !dn/dt where dT/dt=0
  subroutine fcn_tconst(n,tt,x,f)
    use krome_commons
    use krome_ode
    implicit none
    integer::n,ierr
    real*8::x(n),f(n),tt
    call fex(n,tt,x(:),f(:))
    f(idx_Tgas) = 0d0
  end subroutine fcn_tconst

  !*******************************
  subroutine krome_dump(n,rwork,iwork,ni)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_reduction
    use krome_ode
    integer::fnum,i,iwork(:),idx(nrea),j
    real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
    character*16::names(nspec),FMTi,FMTr
    character*50::rnames(nrea),fname,prex
    integer,save::mx_dump=1000 ! max nr of reports before terminating
    fnum = 99
    if (krome_mpi_rank>0) then
      write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
    else
      fname = "KROME_ERROR_REPORT"
    endif
    open(fnum,FILE=trim(fname),status="replace")
    tt = 0d0
    names(:) = get_names()
    rnames(:) = get_rnames()
    call fex(nspec,tt,n(:),dn(:))

    write(fnum,*) "KROME ERROR REPORT"
    write(fnum,*) 
    !SPECIES
    write(fnum,*) "Species abundances"
    write(fnum,*) "**********************"
    write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
    write(fnum,*) "**********************"
    do i=1,nspec
       write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
    end do
    write(fnum,*) "**********************"


    !F90 FRIENDLY RESTART
    write(fnum,*) 
    write(fnum,*) "**********************"
    write(fnum,*) "F90-friendly species"
    write(fnum,*) "**********************"
    do i=1,nspec
       write(prex,'(a,i3,a)') "x(",i,") = "
       write(fnum,*) trim(prex),ni(i),"!"//names(i)
    end do

    write(fnum,*) "**********************"
    
    !RATE COEFFIECIENTS
    k(:) = coe_tab(n(:))
    kmax = maxval(k)
    write(fnum,*)
    write(fnum,*) "Rate coefficients at Tgas",n(idx_Tgas)
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
    write(fnum,*) "**********************"    
    do i=1,nrea
       kperc = 0.d0
       if(kmax>0.d0) kperc = k(i)*1d2/kmax
       write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
    end do
    write(fnum,*) "**********************"
    write(fnum,*)

    !FLUXES
    call load_arrays
    rrmax = fex_check(n(:), n(idx_Tgas))
    idx(:) = idx_sort(arr_flux(:))
    write(fnum,*)
    write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3]"
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
    write(fnum,*) "**********************"    
    do j=1,nrea
       i = idx(j)
       rperc = 0.d0
       if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
       write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
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

    mx_dump = mx_dump - 1
    if (mx_dump==0) stop

  end subroutine krome_dump

  !********************************
  subroutine krome_init()
    use krome_commons
    use krome_tabs
    use krome_subs
    use krome_reduction
#IFKROME_useCoolingZ
    use krome_cooling
#ENDIFKROME
#IFKROME_useStars
    use krome_stars
#ENDIFKROME

    call load_arrays

#IFKROME_useCoolingZ
    call coolingZ_init_tabs()
#ENDIFKROME

#KROME_init_anytab

#IFKROME_useTabs
    call make_ktab()
    call check_tabs()
#ENDIFKROME

#IFKROME_useStars
    call stars_init()
#ENDIFKROME

#IFKROME_useH2esc_omukai
    call init_anytab2D("escape_H2.dat",arrH2esc_ntot(:), &
         arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
         ymulH2esc)
    call test_anytab2D("escape_H2.dat",arrH2esc_ntot(:), &
         arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
         ymulH2esc)
#ENDIFKROME

    !init phys common variables
#KROME_init_phys_variables

    !default for thermo toggle is ON
    !$omp parallel
    krome_thermo_toggle = 1
    !$omp end parallel

    !init photo reactants indexes
#KROME_photopartners

    !get machine precision
    krome_epsilon = epsilon(0d0)

  end subroutine krome_init

  !****************************
  function krome_get_coe(x,Tgas)
    !krome_get_coe: public interface to obtain rate coefficients
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
    real*8::krome_get_coe(nrea),n(nspec),x(:),Tgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
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
