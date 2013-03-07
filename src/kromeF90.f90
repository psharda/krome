module krome_main
contains

#KROME_header

  !********************************
  !KROME main (interface to the solver library)
  subroutine krome(x,rhogas,Tgas,dt)

    use krome_commons
    use krome_subs
    use krome_ode
    USE DVODE_F90_M
    implicit none

    real*8::dt,x(nspec-2),rhogas,Tgas,mass(nspec),n(nspec),tloc,xin
    integer::icount,i

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq,itol,itask,istate,iopt,lrw,liw,mf
#KROME_iwork_array
    real*8::atol(nspec),rtol
#KROME_rwork_array
   
    TYPE (VODE_OPTS) :: OPTIONS

    !****************************
    !init DLSODES (see DLSODES manual)
    !call XSETF(0)!toggle solver verbosity
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0.d0
    itol = 1 !both tolerances are scalar
    rtol = 1d-2 !relative tolerance (default: 1d-4)
    atol = 1d-10 !absolute tolerance (default: 1d-40)
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
    n(1:nspec-2) = x(:)
#ENDIFKROME

    common_ntot = sum(n(1:nspec-2))

    n(idx_Tgas) = Tgas !put temperature in the input array
    nold(:) = n(:) !store initial densities (needed for Jacobian)
    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    do i=1,nspec
       atol(i) = max(n(i)*1d-3, 1d-6)
    end do

    OPTIONS = SET_OPTS(ABSERR_VECTOR=ATOL,RELERR=RTOL,METHOD_FLAG=227)

    do
       icount = icount + 1
       !solve ODE
       !CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
       !     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JES, MF)

       CALL DVODE_F90(FEX, NEQ, n, tloc, dt, ITASK, ISTATE, OPTIONS, &
            J_FCN=JES)

       !check DLSODES exit status
       if(istate==2) then
          exit !sucsessful integration
       elseif(istate==-1) then
          istate = 1 !exceeded internal max iterations
       elseif(istate==-5) then
          istate = 3 !wrong sparsity
       else
          print *,"ERROR: wrong solver exit status!"
          print *,"please check input!"
          print *,"istate:",istate
          stop
       end if
    end do

    !avoid negative species
    do i=1,nspec
       n(i) = max(n(i),0.d0)
    end do

#IFKROME_useX
    x(:) = mass(1:nspec-2)*n(1:nspec-2)/rhogas !return to fractions
    x(:) = x(:) / sum(x) * xin !force mass conservation
#ELSEKROME
    x(:) = n(1:nspec-2)
#ENDIFKROME

    Tgas = n(idx_Tgas) !get new temperature
  end subroutine krome

  !********************************
  subroutine krome_init()
    use krome_tabs
    use krome_subs
#IFKROME_useTabs
    call make_ktab()
    call check_tabs()
#ENDIFKROME
    call load_arrays
  end subroutine krome_init

end module krome_main
