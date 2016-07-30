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
    use krome_getphys
    use krome_ode
    USE DVODE_F90_M
    implicit none

    real*8::dt,x(nmols),rhogas,Tgas,mass(nspec),n(nspec),tloc,xin,ntot
    integer::icount,i

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq,itol,itask,istate,iopt,lrw,liw,mf
#KROME_iwork_array
    real*8::atol(nspec),rtol
#KROME_rwork_array
   
    TYPE (VODE_OPTS) :: OPTIONS
#KROME_iaja_parameters
    integer::iauser(niauser), jauser(njauser)

    !****************************
    !init DLSODES (see DLSODES manual)
    !call XSETF(0)!toggle solver verbosity
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0.d0
    itol = 1 !both tolerances are scalars
    rtol = #KROME_RTOL !relative tolerance
    atol = #KROME_ATOL !absolute tolerance

#KROME_custom_RTOL
#KROME_custom_ATOL

    itask = 1
    iopt = 0

#KROME_maxord

    !MF=
    !  = 227 internal-generated sparse JAC
    !  = 221 user-supplied dense JAC
    !  = 222 internal-generated dense JAC
    !  = 27  user-supplied JAC and sparsity structure
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

    ntot = sum(n(1:nmols))

    n(idx_Tgas) = Tgas !put temperature in the input array

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    OPTIONS = SET_OPTS(ABSERR_VECTOR=ATOL,RELERR=RTOL,METHOD_FLAG=MF,&
         MXHNIL=0)

    CALL USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)

    do
       icount = icount + 1
       !solve ODE
       CALL DVODE_F90(FEX, NEQ, n(:), tloc, dt, ITASK, ISTATE, OPTIONS, &
            J_FCN=JEX)

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
    x(:) = mass(1:nmols)*n(1:nmols)/rhogas !return to fractions
    x(:) = x(:) / sum(x) * xin !force mass conservation
#ELSEKROME
    x(:) = n(1:nmols)
#ENDIFKROME

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !********************************
  subroutine krome_init()
    use krome_tabs
    use krome_subs
    use krome_getphys

    call load_arrays

#IFKROME_useTabs
    call make_ktab()
    call check_tabs()
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
