module krome_main

#IFKROME_useBindC
  use iso_c_binding
#ENDIFKROME_useBindC
  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

#KROME_header

  !*******************************
  !KROME main (interface to the solver library)
#IFKROME_useBindC
#IFKROME_useX
  subroutine krome_c(x,rhogas,Tgas,dt #KROME_dust_arguments #KROME_fexCustom) bind(C,name='krome')
#ELSEKROME
  subroutine krome_c(x,Tgas,dt #KROME_dust_arguments #KROME_fexCustom) bind(C,name='krome')
#ENDIFKROME
    use krome_commons
    use krome_user
    implicit none
    #KROME_double :: Tgas,dt
    #KROME_double :: x(nmols)
#IFKROME_ierr
    integer :: ierr
#ENDIFKROME
#IFKROME_useX
    #KROME_double_value :: rhogas
    call krome(x,rhogas,Tgas,dt #KROME_dust_arguments #KROME_fexCustom)
#ELSEKROME
    real*8 :: rhogas
    call krome(x,Tgas,dt #KROME_dust_arguments #KROME_fexCustom)
#ENDIFKROME
  end subroutine krome_c
#ENDIFKROME_useBindC

#IFKROME_useX
  subroutine krome(x,rhogas,Tgas,dt #KROME_dust_arguments #KROME_fexCustom)
#ELSEKROME
  subroutine krome(x,Tgas,dt #KROME_dust_arguments #KROME_fexCustom)
#ENDIFKROME
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    use krome_getphys
    use krome_tabs
    implicit none
    real*8 :: Tgas,dt
    real*8 :: x(nmols)
#IFKROME_useX
    real*8, value :: rhogas
#ELSEKROME
    real*8 :: rhogas
#ENDIFKROME
    #KROME_externalFexCustom
    real*8::mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
    integer::icount,i,icount_max
    integer:: ierr

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
    rwork(:) = 0d0
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
    n(:) = 0d0 !initialize densities

#IFKROME_useX
    mass(:) = get_mass() !get masses
    xin = sum(x) !store initial fractions
    !compute densities from fractions
    do i = 1,nmols
       if(mass(i)>0d0) n(i) = rhogas * x(i) / mass(i)
    end do
#ELSEKROME
    n(1:nmols) = x(:)
#ENDIFKROME

    n(idx_Tgas) = Tgas !put temperature in the input array

#IFKROME_useDustSizeEvol
    n(nmols+1:nmols+ndust) = krome_dust_asize(:) !set dust sizes
#ENDIFKROME

#IFKROME_usedTdust
    n(nmols+ndust+1:nmols+2*ndust) = krome_dust_T(:)
    call compute_Tdust(n(:),Tgas)
    krome_dust_T(:) = n(nmols+ndust+1:nmols+2*ndust)
#ENDIFKROME

#IFKROME_usePreDustExp
    !pre-calculates exponent
    invTdust(:) = 1d0/(krome_dust_T(:)+1d-40)
    dust_Ebareice_exp(:) = &
         get_Ebareice_exp_array(invTdust(:))
    dust_Ebareice23_exp(:) = &
         get_Ebareice23_exp_array(invTdust(:))
#ENDIFKROME

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

#IFKROME_check_mass_conservation
    mass(:) = get_mass() !get masses
    totmass = sum(n(:) * mass(:)) !calculate total mass
#ENDIFKROME
    !store initial values
    ni(:) = n(:)
    n_global(:) = n(:)

#IFKROME_hasStoreOnceRates
    call makeStoreOnceRates(n(:))
#ENDIFKROME

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
       icount = icount + 1
       !solve ODE
       CALL DLSODES(fex#KROME_postfixFexCustom, NEQ(:), n(:), tloc, dt, &
            ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
            LIW, JES, MF)

#IFKROME_ierr
       !break and return ierr when max count reached
       if(icount>icount_max) then
          ierr = istate
          exit
       end if
#ENDIFKROME

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
       elseif(istate==-3) then
          n(:) = ni(:)
          istate = 1
       else
          got_error = .true.
       end if

       if(got_error.or.icount>icount_max) then
          if (krome_mpi_rank>0) then
            print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
            print *,krome_mpi_rank,"istate:",istate
            print *,krome_mpi_rank,"iter count:",icount
            print *,krome_mpi_rank,"max iter count:",icount_max
            print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
          else
            print *,"ERROR: wrong solver exit status!"
            print *,"istate:",istate
            print *,"iter count:",icount
            print *,"max iter count:",icount_max
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

#KROME_compute_electrons

#IFKROME_check_mass_conservation
    if(abs(1d0-totmass/sum(n(:) * mass(:)))>1d-3) then
       print *,"ERROR: mass conservation failure!"
       print *, tloc,totmass,sum(n(:) * mass(:))
    end if
#ENDIFKROME

    !avoid negative species
    do i=1,nspec
       n(i) = max(n(i),0d0)
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

#IFKROME_useDustSizeEvol
    !returns dust abundance
    krome_dust_asize(:) = n(nmols+1:nmols+ndust)
    krome_dust_asize2(:) = n(nmols+1:nmols+ndust)**2
    krome_dust_asize3(:) = n(nmols+1:nmols+ndust)**3
#ENDIFKROME

#IFKROME_usedTdust
    !returns dust temperature
    krome_dust_T(:) = n(nmols+ndust+1:nmols+2*ndust)
#ENDIFKROME

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !*********************************
  !integrates to equilibrium using constant temperature
#IFKROME_useX
  subroutine krome_equilibrium(x,rhogas,Tgas,verbosity) #KROME_bindC
#ELSEKROME
    subroutine krome_equilibrium(x,Tgas,verbosity) #KROME_bindC
#ENDIFKROME
      use krome_ode
      use krome_subs
      use krome_commons
      use krome_constants
      use krome_getphys
      use krome_tabs
      implicit none
      integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
      integer::i,imax
      integer,optional::verbosity
      integer::verbose
#KROME_double_value :: Tgas
#KROME_double :: x(nmols)
#IFKROME_useX
#KROME_double_value :: rhogas
#ELSEKROME
      real*8 :: rhogas
#ENDIFKROME
      real*8::tloc,n(nspec),mass(nspec),ni(nspec)
      real*8::dt,xin
#KROME_iwork_array
      real*8::atol(nspec),rtol(nspec)
#KROME_rwork_array
      real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
      logical::converged

      integer, save :: ncall=0
      integer, parameter :: ncall_print_frequency=20000
      integer :: ncallp
      integer::charges(nspec)
      real*8::masses(nspec)
      character*16::names(nspec)

      !set verbosity from argument
      verbose = 1 !default is verbose
      if(present(verbosity)) verbose = verbosity

#IFKROME_useCoolingGH
      PLW = 2.11814e-13
      PHI = 1.08928e-13
      PHEI = 2.76947e-14
      PCVI = 1.03070e-17
#ENDIFKROME

      call XSETF(0)!toggle solver verbosity
      meth = 2
      neq = nspec !number of eqns
      liw = size(iwork)
      lrw = size(rwork)
      iwork(:) = 0
      rwork(:) = 0d0
      itol = 4 !both tolerances are scalar
      rtol(:) = 1d-6 !relative tolerance
      atol(:) = 1d-20 !absolute tolerance

      ! Switches to decide when equilibrium has been reached
      ertol = 1d-5  ! relative min change in a species
      eatol = 1d-12 ! absolute min change in a species
      max_time=seconds_per_year*5d8 ! max time we will be integrating for

      !for DLSODES options see its manual
      iopt = 0
      itask = 1
      istate = 1

      mf = 222 !internally evaluated sparsity and jacobian
      tloc = 0d0 !initial time

      n(:) = 0d0 !initialize densities
#IFKROME_useX
      mass(:) = get_mass() !get masses
      xin = sum(x) !store initial fractions
      !compute densities from fractions
      do i = 1,nmols
         if(mass(i)>0d0) n(i) = rhogas * x(i) / mass(i)
      end do
#ELSEKROME
      !copy into array
      n(nmols+1:) = 0d0
      n(1:nmols) = x(:)
#ENDIFKROME

      n(idx_Tgas) = Tgas

      !store previous values
      ni(:) = n(:)
      n_global(:) = ni(:)

#IFKROME_hasStoreOnceRates
      call makeStoreOnceRates(n(:))
#ENDIFKROME

      imax = 1000

      dt = seconds_per_year * 1d2
      t_tot = dt
      converged = .false.
      do while (.not. converged)
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

         !avoid negative species
         do i=1,nspec
            n(i) = max(n(i),0d0)
         end do

#IFKROME_conserve
         n(:) = conserve(n(:),ni(:))
#ENDIFKROME
         ! check if we have converged by comparing the error in any species with an relative abundance above eatol
         converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
              .or. t_tot .gt. max_time

         ! Increase integration time by a reasonable factor
         if(.not. converged) then
            dt = dt * 3.
            t_tot = t_tot + dt
            ni = n
            n_global = n
         endif
      enddo
#IFKROME_useX
      x(:) = mass(1:nmols)*n(1:nmols)/rhogas !return to fractions
#ELSEKROME
      !returns to user array
      x(:) = n(1:nmols)
#ENDIFKROME

      if(t_tot > max_time .and. &
           maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
         print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
         print *, 'Tgas :', Tgas
         names(:) = get_names()
         charges(:) = get_charges()
         masses(:) = get_mass()

         print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
         do i=1,nmols
            print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
         end do
         print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
         print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
         print *, 'Absolute and relative floors:', eatol, ertol
      end if

      ! Print info ever so often
      !$omp critical
      ncall=ncall+1
      ncallp = ncall
      !$omp end critical

      if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
         print *, 'Found equilibrium for ', ncallp, ' cells.'
      end if

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
    use krome_getphys
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
    idx(:) = idx_sort(k(:))
    kmax = maxval(k)
    write(fnum,*)
    write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
    write(fnum,*) "**********************"
    write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
    write(fnum,*) "**********************"
    do j=1,nrea
       i = idx(j)
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
    write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
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
  subroutine krome_init() #KROME_bindC
    use krome_commons
    use krome_tabs
    use krome_subs
    use krome_reduction
    use krome_dust
    use krome_cooling
    use krome_photo
    use krome_fit
    use krome_getphys
#IFKROME_useStars
    use krome_stars
#ENDIFKROME

    !init phys common variables
#KROME_init_phys_variables

    !init metallicity default
    !assuming solar
    total_Z = 1d0

    !default D/D_sol = Z/Z_sol
    !assuming linear scaling
    dust2gas_ratio = total_Z

    !default broadening turubulence velocity
    broadeningVturb2 = 0d0

#IFKROME_useH2dust_constant
    !default clumping factor for
    ! H2 formation on dust by Jura/Gnedin
    clump_factor = 1d0
#ENDIFKROME

    !default for thermo toggle is ON
    !$omp parallel
    krome_thermo_toggle = 1
    !$omp end parallel

    !load arrays with ractants/products indexes
    call load_arrays()

#IFKROME_useH2pd
    !load data for H2 photodissociation (file name default)
    call kpd_H2_loadData()
#ENDIFKROME

#IFKROME_useCoolingZ
    !initialize cooling tabel for metals
    call coolingZ_init_tabs()
#ENDIFKROME

#KROME_init_anytab

#IFKROME_useDustTabs
    call init_dust_tabs()
#ENDIFKROME

#IFKROME_use_GFE_tables
    !init Gibss free energy tables
#KROME_init_GFE_tables

#ENDIFKROME

#IFKROME_useTabs
    call make_ktab()
    call check_tabs()
#ENDIFKROME

#IFKROME_useStars
    call stars_init()
#ENDIFKROME

#IFKROME_useChemisorption
    call init_chemisorption_rates()
#ENDIFKROME

#IFKROME_useH2esc_omukai
    call init_anytab2D("escape_H2.dat",arrH2esc_ntot(:), &
         arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
         ymulH2esc)
    call test_anytab2D("escape_H2.dat",arrH2esc_ntot(:), &
         arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
         ymulH2esc)
#ENDIFKROME

#IFKROME_useMayerOpacity
    !call init_anytab2D("mayer_E2.dat",mayer_x(:), &
    !     mayer_y(:), mayer_z(:,:), mayer_xmul, &
    !     mayer_ymul)
    !call test_anytab2D("mayer_E2.dat",mayer_x(:), &
    !     mayer_y(:), mayer_z(:,:), mayer_xmul, &
    !     mayer_ymul)
#ENDIFKROME

#IFKROME_useCoolingCO
    !initialize CO cooling
    call init_coolingCO()
#ENDIFKROME

#IFKROME_useCoolingOH
    !initialize OH cooling
    call init_coolingOH()
#ENDIFKROME

#IFKROME_useCoolingH2O
    !initialize H2O cooling
    call init_coolingH2O()
#ENDIFKROME

#IFKROME_useCoolingHCN
    !initialize HCN cooling
    call init_coolingHCN()
#ENDIFKROME

#IFKROME_useCoolingZCIE
    !initialize metal CIE cooling
    call init_coolingZCIE()
#ENDIFKROME

#IFKROME_useCoolingZCIENOUV
  !initialize metal CIE cooling no UV case
  call init_anytab2D("coolZ_CIE2012NOUV.dat",CoolZNOUV_x(:), &
        CoolZNOUV_y(:), CoolZNOUV_z(:,:), CoolZNOUV_xmul, &
        CoolZNOUV_ymul)
  call test_anytab2D("coolZ_CIE2012NOUV.dat",CoolZNOUV_x(:), &
        CoolZNOUV_y(:), CoolZNOUV_z(:,:), CoolZNOUV_xmul, &
        CoolZNOUV_ymul)
#ENDIFKROME

    !initialize the table for exp(-a/T) function
    call init_exp_table()

#IFKROME_useGammaPop
    call load_parts()
#ENDIFKROME

    !init photo reactants indexes
#KROME_photopartners

    !get machine precision
    krome_epsilon = epsilon(0d0)

#IFKROME_useSemenov
    !get binding energies 
    Ebinding(:) = get_EbindBare()
#ENDIFKROME
    !load verbatim reactions
    call loadReactionsVerbatim()

  end subroutine krome_init

  !****************************
  function krome_get_coe(x,Tgas) #KROME_bindC
    !krome_get_coe: public interface to obtain rate coefficients
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
#IFKROME_useBindC
    real(kind=c_double) :: x(nmols)
    real(kind=c_double), value :: Tgas
    real(kind=c_double), target :: coes(nrea)
    type(c_ptr) :: krome_get_coe
#ELSEKROME_useBindC
    real*8 :: krome_get_coe(nrea), x(nmols), Tgas
#ENDIFKROME_useBindC
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
#IFKROME_useBindC
    coes(:) = coe_tab(n(:))
    krome_get_coe = c_loc(coes)
#ELSEKROME_useBindC
    krome_get_coe(:) = coe_tab(n(:))
#ENDIFKROME_useBindC

  end function krome_get_coe

  !****************************
  function krome_get_coeT(Tgas) #KROME_bindC
    !krome_get_coeT: public interface to obtain rate coefficients
    ! with argument Tgas only
    use krome_commons
    use krome_subs
    use krome_tabs
    implicit none
#IFKROME_useBindC
    real(kind=c_double), value :: Tgas
    real(kind=c_double), target :: coeTs(nrea)
    type(c_ptr) :: krome_get_coeT
#ELSEKROME_useBindC
    real*8 :: krome_get_coeT(nrea),Tgas
#ENDIFKROME_useBindC
    real*8::n(nspec)
    n(idx_Tgas) = Tgas
#IFKROME_useBindC
    coeTs(:) = coe_tab(n(:))
    krome_get_coeT = c_loc(coeTs)
#ELSEKROME_useBindC
    krome_get_coeT(:) = coe_tab(n(:))
#ENDIFKROME_useBindC
  end function krome_get_coeT


#IFKROME_reducer
  !****************************
  !this subroutine randomly performs integrations over the
  ! intervals of the variables provided via the arguments.
  ! The runs are employed to track the most active reactions
  subroutine krome_reducer(#KROME_reducerVarsInterface,&
       verbosity,iterations,treshold)
    use krome_commons
    use krome_constants
    use krome_user
    use krome_subs
#KROME_useIFPORT
    implicit none
    integer::imax,i,verbose,idx(nrea),j
    integer,optional::verbosity,iterations
    character*50::name(nrea),fname
    real*8,optional::treshold
    real*8::dt,t,tmax,x(nmols),flux(nrea),fluxmax,frac,points(nrea)
#KROME_reducerVarsDeclare

    fname = "KROME_reduced.dat"

    !option variable for verbosity
    if(present(verbosity)) then
       verbose = verbosity
    else
       verbose = 2
    end if

    !option variable for iterations
    if(present(iterations)) then
       imax = iterations
    else
       imax = 10
    end if

    !option variable for treshold
    if(present(treshold)) then
       frac = treshold
    else
       frac = 1d-6
    end if

    if(frac.ge.1d0) then
       print *,"ERROR: the fraction in reducer is >=1!"
       stop
    end if

    if(verbose>0) then
       print *,"reducer starts!"
       print *,"threshold fraction:",frac
       print *,"iterations:",imax
    end if

    !init points to zero
    points(:) = 0d0

#KROME_reducerVarsLog

    tmax = 1d7*seconds_per_year !max time (s)
    do i=1,imax

#KROME_reducerVarsRandomize

#KROME_reducerVarsUserSet

       if(verbose>1) print *,"***************"
       if(verbose>0) print *,i,"of",imax
       if(verbose>1) then
#KROME_reducerPrintInits
       end if

       !initialize variables
       x(:) = 1d-40
       x(idx_H) = ntot
       call krome_scale_Z(x(:),Zmetals)
       t = 0d0
       dt = 1d-3*seconds_per_year
       do
          dt = dt * 1.1
          t = t + dt
          call krome(x(:),Tgas,dt)
          flux(:) = krome_get_flux(x(:),Tgas) !get fluxes
          fluxmax = maxval(flux)
          !call the sorting algorithm (bubblesort)
          do j=1,nrea
             if(flux(j)>fluxmax*frac) points(j) = points(j) + 1d0
          end do
          if(t>tmax) exit
       end do
    end do

    idx(:) = idx_sort(points(:))
    name(:) = get_rnames() !get reaction names

    open(32,file=trim(fname),status="replace")
    write(32,*) "**** RESULTS ****"
    write(32,'(a50,a8,a10)') "reaction","idx","score"
    do j=1,size(idx)
       write(32,'(a50,I8,F10.0)') name(idx(j)),idx(j),points(idx(j))
    end do
    close(32)
    if(verbose>0) then
       print *,"reducer ends!"
       print *,"data saved to "//trim(fname)
    end if
  end subroutine krome_reducer

#ENDIFKROME

end module krome_main
