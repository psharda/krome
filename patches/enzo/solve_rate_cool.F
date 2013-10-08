#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, fh, dtoh, z_solar,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                ndratec, dtemstart, dtemend, h2dusta,
     &                ncrna, ncrd1a, ncrd2a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela, 
     &                gasgra, metala, n_xe, xe_start, xe_end,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, iradcoupled, iradstep, ierr,
     &                irt_honly,
     &                kphHI, kphHeI, kphHeII, kdissH2I, 
     &                photogamma, 
     &                ih2optical, iciecool, ithreebody, ciecoa, 
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating)

!
!  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
!  modified2:  October, 1996 by GB; adapted to AMR
!  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
!  modified4:  June,    2005 by GB to solve rate & cool at same time
!  modified5:  April,   2009 by JHW to include radiative transfer
!  modified6:  September, 2009 by BDS to include cloudy cooling
!  modified7:  2013 by KROME developers to include KROME package
!
!  PURPOSE:
!    Solve the multi-species rate and cool equations.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!    HI,HII   - H density fields (neutral & ionized)
!    HeI/II/III - He density fields
!    DI/II    - D density fields (neutral & ionized)
!    HDI      - neutral HD molecule density field
!    HM       - H- density field
!    H2I      - H_2 (molecular H) density field
!    H2II     - H_2+ density field
!    kph*     - photoionization fields
!    gamma*   - photoheating fields
!
!    is,ie    - start and end indices of active region (zero based)
!    idual    - dual energy formalism flag (0 = off, 1 = on)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    idim     - dimensionality (rank) of problem
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    iradshield - flag for crude radiative shielding correction 
!    iradtype - type of radiative field (only used if = 8)
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!    imcool   - flag if there is metal cooling
!    idust    - flag for H2 formation on dust grains
!    imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
!    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
!    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
!    iradtrans - flag to include radiative transfer (1 = on, 0 = off)
!    iradcoupled - flag to indicate coupled radiative transfer
!    iradstep - flag to indicate intermediate coupled radiative transfer timestep
!
!    fh       - Hydrogen mass fraction (typically 0.76)
!    dtoh     - Deuterium to H mass ratio
!    z_solar  - Solar metal mass fraction
!    dt       - timestep to integrate over
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!    nratec   - dimensions of chemical rate arrays (functions of temperature)
!    dtemstart, dtemend - start and end of dust temperature range
!    ndratec  - extra dimension for H2 formation on dust rate (dust temperature)
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3, clPar4, clPar5 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!
!  OUTPUTS:
!    update chemical rate densities (HI, HII, etc)
!
!  PARAMETERS:
!    itmax   - maximum allowed sub-cycle iterations
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------
!     USE KROME
      use krome_main
      implicit NONE
#include "fortran_types.def"

!  General Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, iradtrans,
     &        iradcoupled, iradstep, n_xe, ierr, imcool, idust,
     &        irt_honly, igammah, ih2optical, iciecool, ithreebody,
     &        ndratec
      R_PREC    dt, aye, temstart, temend, eta1, eta2, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh, z_solar, 
     &        xe_start, xe_end, dtemstart, dtemend

!  Density, energy and velocity fields fields

      R_PREC    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      R_PREC    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn)

!  Radiation fields

      R_PREC kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &     kdissH2I(in,jn,kn)
      R_PREC photogamma(in,jn,kn)

!  Cooling tables (coolings rates as a function of temperature)

      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      R_PREC    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec), 
     &        ciecoa(nratec)
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    metala(nratec, n_xe)
      R_PREC    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
      R_PREC    gammaha 

!  Chemistry tables (rates as a function of temperature)

      R_PREC k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &     k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7), h2dusta(nratec, ndratec),
     &     ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &     k24, k25, k26, k27, k28, k29, k30, k31

!  Cloudy cooling data

      INTG_PREC icmbTfloor, iClHeat, clGridRank, clDataSize
      INTG_PREC clGridDim(5)
      R_PREC clEleFra
      R_PREC clPar1(clGridDim(1)), clPar2(clGridDim(2)), 
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)), 
     &     clPar5(clGridDim(5))
      R_PREC clCooling(clDataSize), clHeating(clDataSize)

!  Parameters

      INTG_PREC itmax, ijk
      parameter (itmax = 10000, ijk = MAX_ANY_SINGLE_DIRECTION)

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

!  Locals

      INTG_PREC i, j, k, iter
      INTG_PREC clGridDim1, clGridDim2, clGridDim3, clGridDim4, 
     &     clGridDim5
      R_PREC ttmin, dom, energy, comp1, comp2
      real*8 coolunit, dbase1, tbase1, xbase1, chunit, uvel
      real*8 heq1, heq2, eqk221, eqk222, eqk131, eqk132,
     &                 eqt1, eqt2, eqtdef, dheq, heq, dlogtem

!  row temporaries

      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), 
     &     dtit(ijk), ttot(ijk), p2d(ijk), tgas(ijk), tgasold(ijk),
     &     tdust(ijk), metallicity(ijk), rhoH(ijk), olddtit

!  Rate equation row temporaries

      R_PREC HIp(ijk), HIIp(ijk), HeIp(ijk), HeIIp(ijk), HeIIIp(ijk),
     &     HMp(ijk), H2Ip(ijk), H2IIp(ijk),
     &     dep(ijk), dedot(ijk),HIdot(ijk), dedot_prev(ijk),
     &     DIp(ijk), DIIp(ijk), HDIp(ijk), HIdot_prev(ijk),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)
      R_PREC k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k50(ijk), k51(ijk), k52(ijk), k53(ijk), k54(ijk),
     &     k55(ijk), k56(ijk), k13dd(ijk, 7), h2dust(ijk),
     &     ncrn(ijk), ncrd1(ijk), ncrd2(ijk)

!  Cooling/heating row locals

      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk), 
     &     brem(ijk), edot(ijk)
      R_PREC hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk), 
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk), cieco(ijk)

      R_PREC nHI, nH2I, atten, H2delta

!  Iteration mask

      LOGIC_PREC itmask(ijk)

!  New variables KROME
      !NSPECIES IS IN FORTRAN.DEF
      INTG_PREC nsp
      parameter (nsp=9)
      real*8 ugox(nsp)
      real*8 rhogas
      real*8 dt_hydro
      real*8 time_conversion
      parameter (time_conversion=3.1557e+13)

!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set error indicator

      ierr = 0
      
!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      uvel     = uxyz / utim
c      chunit = (7.17775e-12_RKIND)/(2._RKIND*uvel*uvel*mh)   ! 4.5 eV per H2 formed
      chunit = (1.60218e-12_RKIND)/(2._RKIND*uvel*uvel*mh)   ! 1 eV per H2 formed

      dlogtem = (log(temend) - log(temstart))/REAL(nratec-1,RKIND)

!  Convert densities from comoving to proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**(-3))

      call ceiling_species(d, de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal)

!  Loop over zones, and do an entire i-column in one go

      do k = ks+1, ke+1
      do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

!       Set iteration mask to include only cells with radiation in the
!       intermediate coupled chemistry / energy step

         if (iradcoupled == 1 .and. iradstep == 1) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) > 0) then 
                  itmask(i) = .true. 
               else 
                  itmask(i) = .false.
               endif
            enddo
         endif

!     Normal rate solver, but don't double-count the cells with radiation
 
         if (iradcoupled == 1 .and. iradstep == 0) then
            do i = is+1, ie+1
               if (kphHI(i,j,k) > 0) then 
                  itmask(i) = .false. 
               else 
                  itmask(i) = .true.
               endif
            enddo
         endif

!     No radiation timestep coupling

         if (iradcoupled == 0 .or. iradtrans == 0) then 
            do i = is+1, ie+1
               itmask(i) = .true.
            enddo
         endif

!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0._RKIND
         enddo

         do i = is+1, ie+1

            !NEW ROUTINE TO EVALUATE THE TEMPERATURE FROM INTERNAL ENERGY
            !KROME 2013
            call evaluate_temp(d,e,ge,u,
     &                v,w,de,HI,HII,
     &                HeI,HeII,HeIII,
     &                idual,imethod,in, jn, kn,is, ie, j, k,
     &                iexpand,ispecies,imetal,idim, 
     &                aye,temstart,temend,z_solar,
     &                utem,uxyz,uaye,urho,utim,gamma,
     &                HM,H2I,H2II,
     &                DI,DII,HDI,metal,
     &                tgas,tgasold,p2d,metallicity,rhoH,
     &                itmask)

              !NEW SPECIES ARRAY IN NUMBER DENSITY
              !TO UPDATE IF USE NEW SPECIES KROME 2013
              ugox(1)=HI(i,j,k)*dom
              ugox(2)=de(i,j,k)*dom
              ugox(3)=HII(i,j,k)*dom
              ugox(4)=HeI(i,j,k)*dom/4.0d0 
              ugox(5)=HeII(i,j,k)*dom/4.0d0
              ugox(6)=HeIII(i,j,k)*dom/4.0d0
              ugox(7)=HM(i,j,k)*dom
              ugox(8)=H2I(i,j,k)*dom/2.0d0
              ugox(9)=H2II(i,j,k)*dom/2.0d0
              !ugox(11)=DI(i,j,k)*dom 
              !ugox(12)=DII(i,j,k)*dom
              !ugox(13)=HDI(i,j,k)*dom/3.0d0
              
              !CONVERT FROM CODE UNITS TO SECONDS AND 1/cm3 KROME 2013
              rhogas=d(i,j,k)*dom
              dt_hydro= utim*dt !dt*time_conversion

              !CALL SOLVER: COOLING/HEATING COUPLED WITH SOLVER KROME 2013
              call krome(ugox(:),rhogas,tgas(i),dt_hydro) 

              !CONVERT ABUNDANCES FROM NUMBER DENSITY TO CODE UNITS KROME 2013
              HI(i,j,k)   = ugox(1)/dom
              de(i,j,k)   = ugox(2)/dom
              HII(i,j,k)  = ugox(3)/dom
              HeI(i,j,k)  = ugox(4)*4.0d0/dom 
              HeII(i,j,k) = ugox(5)*4.0d0/dom
              HeIII(i,j,k)= ugox(6)*4.0d0/dom
              HM(i,j,k)   = ugox(7)/dom
              H2I(i,j,k)  = ugox(8)*2.0d0/dom
              H2II(i,j,k) = ugox(9)*2.0d0/dom
              !DI(i,j,k)   = ugox(11)/dom 
              !DII(i,j,k)  = ugox(12)/dom
              !HDI(i,j,k)  = ugox(13)/dom*3.0d0

              !EVALUATE ENERGY VARIATION FROM THE OLD
              !AND NEW TEMPERATURES KROME 2013
              edot(i)=(tgas(i)-tgasold(i))*d(i,j,k)
     &             /((gamma-1.0d0)*utem*dt)
           
            enddo  

            !UPDATE INTERNAL ENERGY AND TOTAL ENERGY KROME 2013
            do i = is+1, ie+1
               if (itmask(i)) then
                 e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dt
                 if (idual .eq. 1) then
                    ge(i,j,k) = ge(i,j,k)+ edot(i)/d(i,j,k)*dt
                endif
              endif
            enddo

!     Next j,k

       enddo
      enddo

!     Convert densities back to comoving from proper

      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**3)

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, 
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh)

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
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
     &                           HM, H2I, H2II, DI, DII, HDI, metal,
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
      R_PREC metal(in,jn,kn)

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
