#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!/////////////////////  SUBROUTINE KROME_DRIVER  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine krome_driver(d, e, ge, u, v, w, de,
#KROME_driver_args
     &     metal,in, jn, kn, imethod,
     &     idual, ispecies, imetal, idim,
     &     is, js, ks, ie, je, ke,
     &     dt, aye, temstart, 
     &     utem, uxyz, uaye, urho, utim,
     &     gamma, fh, dtoh, z_solar, ierr)

!     
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

      implicit NONE
#include "fortran_types.def"

!     General Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, imethod,
     &     idual, ispecies, imetal, idim,
     &     ierr
      R_PREC    dt, aye, temstart, gamma,
     &     utim, uxyz, uaye, urho, utem, fh, dtoh, z_solar

!     Density, energy and velocity fields fields
#KROME_rprec

      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &     u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &     metal(in,jn,kn)


!     Parameters
      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.e-5_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-10_RKIND)
#endif

      real*8 mh
      parameter (mh = mass_h)

!     Locals

      INTG_PREC i, j, k
      R_PREC dom

!     row temporaries

      R_PREC tgas(ijk), tgasold(ijk), rhoH(ijk), metallicity(ijk)

!     Rate equation row temporaries

      R_PREC edot(ijk)

!     New variables for KROME
      INTG_PREC nsp
      parameter (nsp=10)
      real*8 krome_x(nsp)
      real*8 rhogas
      real*8 dt_hydro,factor
      real*8 time_conversion
      parameter (time_conversion=3.1557e+13)


!     
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set error indicator

      ierr = 0
      
!     Set units

      dom      = urho*(aye**3)/mh

!     Convert densities from comoving to proper
!     ORIGINAL ENZO ROUTINE, IF YOU ADDITIONAL SPECIES?!
      factor = aye**(-3) !scaling factor for comoving->proper

!     Loop over zones, and do an entire i-column in one go
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
!     scale comoving->proper
               d(i,j,k) = d(i,j,k)*factor
#KROME_scale

!     check minimal value
#KROME_minval

!     HERE WE SHOULD CALL THE EVALUATE TEMPERATURE ROUTINE
!     THIS CAN BE ALSO RE-WRITTEN
               call evaluate_temp(
     &              d,e,ge,u,v,w,de,HI,HII,HeI,HeII,HeIII,
     &              HM,H2I,H2II,DI,DII,HDI,HeHI,metal,
     &              idual,imethod,in, jn, kn,is, ie, j, k,
     &              ispecies,imetal,idim, 
     &              aye,temstart,z_solar,
     &              utem,uxyz,uaye,urho,utim,gamma,
     &              tgas,tgasold,metallicity,rhoH)

!NEW SPECIES ARRAY IN NUMBER DENSITY
!TO UPDATE IF USE NEW SPECIES KROME 2013
               krome_x(1)=HI(i,j,k)*dom
               krome_x(2)=de(i,j,k)*dom
               krome_x(3)=HII(i,j,k)*dom
               krome_x(4)=HeI(i,j,k)*dom/4.0d0 
               krome_x(5)=HeII(i,j,k)*dom/4.0d0
               krome_x(6)=HeIII(i,j,k)*dom/4.0d0
               krome_x(7)=HM(i,j,k)*dom
               krome_x(8)=H2I(i,j,k)*dom/2.0d0
               krome_x(9)=H2II(i,j,k)*dom/2.0d0
               krome_x(10)=HeHI(i,j,k)*dom/5.0d0 !TO BE CHECKED
!krome_x(11)=DI(i,j,k)*dom 
!krome_x(12)=DII(i,j,k)*dom
!krome_x(13)=HDI(i,j,k)*dom/3.0d0
               
!CONVERT FROM CODE UNITS TO SECONDS AND 1/cm3 KROME 2013
               rhogas=d(i,j,k)*dom
               dt_hydro= utim*dt !dt*time_conversion

!CALL SOLVER: COOLING/HEATING COUPLED WITH SOLVER KROME 2013
               call krome(krome_x(:),tgas(i),dt_hydro) 

!CONVERT ABUNDANCES FROM NUMBER DENSITY TO CODE UNITS KROME 2013
               HI(i,j,k)   = krome_x(1)/dom
               de(i,j,k)   = krome_x(2)/dom
               HII(i,j,k)  = krome_x(3)/dom
               HeI(i,j,k)  = krome_x(4)*4.0d0/dom 
               HeII(i,j,k) = krome_x(5)*4.0d0/dom
               HeIII(i,j,k)= krome_x(6)*4.0d0/dom
               HM(i,j,k)   = krome_x(7)/dom
               H2I(i,j,k)  = krome_x(8)*2.0d0/dom
               H2II(i,j,k) = krome_x(9)*2.0d0/dom
               HeHI(i,j,k) = krome_x(10)*5.0d0/dom !TO BE CHECKED
!DI(i,j,k)   = krome_x(11)/dom 
!DII(i,j,k)  = krome_x(12)/dom
!HDI(i,j,k)  = krome_x(13)/dom*3.0d0

!EVALUATE ENERGY VARIATION FROM THE OLD
!AND NEW TEMPERATURES KROME 2013
!edot = de/dt
!dT / de = utem * (gamma - 1.0) / d
!de/dt = dT/dt * d / (gamma - 1.0) / utem
               edot(i)=(tgas(i)-tgasold(i))*d(i,j,k)
     &              /((gamma-1.0d0)*utem*dt)
               
            enddo

!UPDATE INTERNAL ENERGY AND TOTAL ENERGY KROME 2013
            do i = is+1, ie+1
               e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dt
               if (idual .eq. 1) then
                  ge(i,j,k) = ge(i,j,k)+ edot(i)/d(i,j,k)*dt
               endif
            enddo


!     Next j,k
!     
         enddo
      enddo

!     Convert densities back to comoving from proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI, HeHI, metal,
     &     is, ie, js, je, ks, ke,
     &     in, jn, kn, ispecies, imetal, aye**3)

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI, metal, 
     &     d, is, ie, js, je, ks, ke,
     &     in, jn, kn, ispecies, imetal, fh, dtoh)

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, HeHI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC    HeHI(in,jn,kn)
      R_PREC     d(in,jn,kn),metal(in,jn,kn)
      R_PREC    factor

!     locals

      INTG_PREC i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
               de(i,j,k)    = de(i,j,k)*factor
               HI(i,j,k)    = HI(i,j,k)*factor
               HII(i,j,k)   = HII(i,j,k)*factor
               HeI(i,j,k)   = HeI(i,j,k)*factor
               HeII(i,j,k)  = HeII(i,j,k)*factor
               HeIII(i,j,k) = HeIII(i,j,k)*factor
               HeHI(i,j,k)  = HeHI(i,j,k)*factor !SB2013
            enddo
         enddo
      enddo
      if (ispecies > 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies > 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (imetal == 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
               enddo
            enddo
         enddo
      endif

      return
      end

c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species(d, de, HI, HII, HeI, HeII, HeIII,
     &                           HM, H2I, H2II, DI, DII, HDI,HeHI,metal,
     &                           is, ie, js, je, ks, ke,
     &                           in, jn, kn, ispecies, imetal)
c -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC     d(in,jn,kn)
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC metal(in,jn,kn), HeHI(in,jn,kn) !SB2013

!     locals

      INTG_PREC i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               de(i,j,k)    = max(de(i,j,k), tiny)
               HI(i,j,k)    = max(HI(i,j,k), tiny)
               HII(i,j,k)   = max(HII(i,j,k), tiny)
               HeI(i,j,k)   = max(HeI(i,j,k), tiny)
               HeII(i,j,k)  = max(HeII(i,j,k), tiny)
               HeIII(i,j,k) = max(HeIII(i,j,k), 1.e-5_RKIND*tiny)
               HeHI(i,j,k)  = max(HeHI(i,j,k), tiny)
            enddo
         enddo
      enddo
      if (ispecies > 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = max(HM(i,j,k), tiny)
                  H2I(i,j,k)  = max(H2I(i,j,k), tiny)
                  H2II(i,j,k) = max(H2II(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies > 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = max(DI(i,j,k), tiny)
                  DII(i,j,k) = max(DII(i,j,k), tiny)
                  HDI(i,j,k) = max(HDI(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (imetal == 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = min(max(metal(i,j,k),
     &                               tiny), 0.9_RKIND*d(i,j,k))
               enddo
            enddo
         enddo
      endif

      return
      end


! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, fh, dtoh)
! -------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!     Arguments

      INTG_PREC in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn), metal(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      R_PREC    fh, dtoh

!     Parameters

      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)

!     locals

      INTG_PREC i, j, k
      R_PREC totalH(ijk), totalHe(ijk), 
     &     totalD, correctH, correctHe, correctD
      R_PREC metalfree(ijk)

!     Loop over all zones

      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      if (imetal == 1) then
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k) - metal(i,j,k)
         enddo
      else
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k)
         enddo
      endif

      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
      enddo

!     include molecular hydrogen

      if (ispecies > 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
         enddo
      endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
         correctH = fh*metalfree(i)/totalH(i)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH

         correctHe = (1._RKIND - fh)*metalfree(i)/totalHe(i)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies > 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies > 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + 
     &           2._RKIND/3._RKIND*HDI(i,j,k)
            correctD = fh*dtoh*metalfree(i)/totalD
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND + 
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies > 1) de(i,j,k) = de(i,j,k)
     &             - HM(i,j,k) + H2II(i,j,k)/2._RKIND
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end
