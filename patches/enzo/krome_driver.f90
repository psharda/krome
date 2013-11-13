  #include "fortran.def"
  #include "phys_const.def"
  #include "error.def"
  
  !KROME_DRIVER
  
  subroutine krome_driver(d, e, ge, u, v, w, &
       #KROME_args
       in, jn, kn, imethod, &
       idual, idim, &
       is, js, ks, ie, je, ke, &
       dt, aye, temstart, &
       utem, uxyz, uaye, urho, utim, & 
       gamma, fh, dtoh)
    
    
    !     SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
    !     RE-WRITTEN FROM THE ORIGINAL ENZO SOLVE_RATE SUBROUTINE 
    !     
    !     2013, KROME DEVELOPERS to interface the package with ENZO
    !     
    !     PURPOSE:
    !     Solve the multi-species rate and cool equations via KROME.
    !     
    !     INPUTS:
    !     in,jn,kn - dimensions of 3D fields
    !     
    !     d        - total density field
    !     de       - electron density field
    !     HI,HII   - H density fields (neutral & ionized)
    !     HeI/II/III - He density fields
    !     DI/II    - D density fields (neutral & ionized)
    !     HDI      - neutral HD molecule density field
    !     HM       - H- density field
    !     H2I      - H_2 (molecular H) density field
    !     H2II     - H_2+ density field
    !     
    !     is,ie    - start and end indices of active region (zero based)
    !     idual    - dual energy formalism flag (0 = off, 1 = on)
    !     idim     - dimensionality (rank) of problem
    !     ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
    !     imetal   - flag if metal field is active (0 = no, 1 = yes)
    !     imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
    !     temstart - start of temperature range for rate table
    !     
    !     fh       - Hydrogen mass fraction (typically 0.76)
    !     dtoh     - Deuterium to H mass ratio
    !     z_solar  - Solar metal mass fraction
    !     dt       - timestep to integrate over
    !     aye      - expansion factor (in code units)
    !     
    !     utim     - time units (i.e. code units to CGS conversion factor)
    !     uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
    !     urho     - density units
    !     uxyz     - length units
    !     utem     - temperature(-like) units
    !     
    !     OUTPUTS:
    !     update chemical rate densities (HI, HII, etc) and energy
    !     
    !     PARAMETERS:
    !     mh      - H mass in cgs units
    !     
    !-----------------------------------------------------------------------
    !     USE KROME
    use krome_main
    use krome_user

    implicit none
#include "fortran_types.def"

    real*8,parameter::mh=mass_h
    real*8::dom,,factor,tgas,tgasold
    real*8::utim,dt,dt_hydro,rhogas,idom,edot,krome_tiny
    R_PREC d(:,:,:),e(:,:,:),ge(:,:,:)
    R_PREC u(:,:,:),v(:,:,:),w(:,:,:)
    R_PREC dt,aye,temstart,utem,uxyz,uaye,urho,utim,gamma
    INTG_PREC in,jn,kn,imethod,idual,is,js,ks,ie,je,ke
    integer::i,j,k

#KROME_rprec

    !******************************

    !set error indicator
    ierr = 0

    !set units
    dom = urho*(aye**3)/mh

    !scaling factor for comoving->proper
    factor = aye**(-3)

    !scale comoving->proper
    d(:,:,:) = d(:,:,:) * factor
#KROME_scale

    !check minimal value
    krome_tiny = tiny
    do k = ks+1, ke+1
       do j = js+1, je+1
          do i = is+1, ie+1  
#KROME_minval
          end do
       end do
    end do

    !loop over zones
    do k = ks+1, ke+1
       do j = js+1, je+1
          do i = is+1, ie+1

             rhogas = #KROME_sum

             !convert to number densities
#KROME_dom
             call evaluate_temp(d(i,j,k), e(i,j,k), ge(i,j,k),&
                  u(i,j,k), v(i,j,k), w(i,j,k),& 
                  krome_x(:),imethod,idual,idim,tgas,&
                  temstart,utem,rhogas)

             !store old tgas
             tgasold = tgas

             !convert to g/cm3
             rhogas = d(i,j,k) * dom
             dt_hydro = utim*dt !dt*time_conversion

             !call KROME solver
             call krome(krome_x(:),tgas,dt_hydro) 

             idom = 1.d0/dom
             !convert back to code units
#KROME_mod

             !evaluate energy from temperature difference
             edot = (tgas - tgasold) * d(i,j,k) &
                  / ((gamma - 1.d0) * utem * dt)

             !update internal energy
             e(i,j,k)  = e(i,j,k) + edot / d(i,j,k) * dt
             !when using dual
             if (idual .eq. 1) ge(i,j,k) = ge(i,j,k)+ edot / d(i,j,k) * dt
          end do
       end do
    end do

    !scale comoving<-proper
    factor = aye**3
    d(:,:,:) = d(:,:,:) * factor
#KROME_scale

  end subroutine krome_driver
