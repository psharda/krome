#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!//////////////////////  SUBROUTINE EVALUATE_TEMP  \\\\\\\\\\\\\\\\\\\\\\

      subroutine evaluate_temp(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                idual, imethod,in, jn, kn,is, ie, j, k,
     &                iexpand, ispecies, imetal, idim, 
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim, gamma,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                tgas, tgasold, p2d, metallicity, rhoH, 
     &                itmask)

!  EVALUATE TEMPERATURE FROM ENERGY FOR A GIVEN CELL 
!  written by: KROME DEVELOPERS 
!  date: November, 2012
!
!  PURPOSE:
!    Evaluate temperature Tgas from the internal energy e(i,j,k)
!    or ge(i,j,k) depending on the energy flag idual
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!  Arguments

      INTG_PREC in, jn, kn, is, ie, j, k
      INTG_PREC imethod, idim,
     &        idual, iexpand, ispecies, imetal
      R_PREC    aye, temstart, temend, z_solar,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma

      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      R_PREC tgas(in), tgasold(in)

!  Parameters

      real*8 mh
      R_PREC    ZSOLAR, mu_metal
      parameter (mh = mass_h)         !DPC
      parameter (mu_metal = 16.d0)    ! approx. mean molecular weight of metals


!  Locals

      INTG_PREC i, j1, iradfield
      R_PREC dom, energy
      R_PREC gamma2, x, dom_inv
      R_PREC coolunit, dbase1, tbase1, xbase1,
     &     nH2, nother

!  Slice locals
 
      R_PREC p2d(in), rhoH(in), nh(in), metallicity(in)
!  Metal cooling locals

      R_PREC xe, log_xe, xe1, dlogT, xi
      R_PREC metalfree

!  Iteration mask

      LOGIC_PREC itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set units

      dom      = urho*(aye**3)/mh
      dom_inv  = 1.d0/dom
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      iradfield = -1

!     Set compton cooling coefficients (and temperature)
!
!      if (iexpand .eq. 1) then
!         comp1 = compa  * (1.d0 + zr)**4
!         comp2 = 2.73d0 * (1.d0 + zr)
!      else
!         comp1 = 1d-35  ! tiny was too large in comparison with compa
!         comp2 = tiny
!      endif

!     Compute Pressure

      if (imethod .eq. 2) then

!        Zeus - e() is really gas energy

         do i = is+1, ie+1
          if(itmask(i)) then
             p2d(i) = (gamma - 1.d0)*d(i,j,k)*e(i,j,k)
          end if
         enddo
      else
          if(idual.eq.1) then

!           PPM with dual energy -- use gas energy

            do i = is+1, ie+1
             if(itmask(i)) then
               p2d(i) = (gamma - 1.d0)*d(i,j,k)*ge(i,j,k)
             end if
            enddo
          else

!           PPM without dual energy -- use total energy

            do i = is+1, ie+1
             if(itmask(i)) then
               p2d(i) = e(i,j,k) - 0.5d0*u(i,j,k)**2
               if (idim.gt.1) p2d(i) = p2d(i) - 0.5d0*v(i,j,k)**2
               if (idim.gt.2) p2d(i) = p2d(i) - 0.5d0*w(i,j,k)**2
               p2d(i) = max((gamma - 1.d0)*d(i,j,k)*p2d(i), tiny)
             end if
            enddo
          endif
      endif

!     Compute temperature and H number density

      do i = is+1, ie+1
         if (itmask(i)) then
            tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            rhoH(i) = HI(i,j,k) + HII(i,j,k)
         end if
      enddo

!          (include molecular hydrogen, but ignore deuterium)

      if (ispecies.gt.1) then
         do i = is+1, ie+1
         if (itmask(i)) then
             tgas(i) = tgas(i) +
     &            HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.d0
             rhoH(i) = rhoH(i) + H2I(i,j,k) + H2II(i,j,k)
          end if
          enddo
      endif

!     Include metal species and calculate metallicity
      
      if (imetal.eq.1) then
         do i = is+1, ie+1
         tgas(i) = tgas(i) + metal(i,j,k)/mu_metal
         metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
         enddo
      endif

      do i = is+1, ie+1
      if (itmask(i)) then
         tgas(i) = max(p2d(i)*utem/tgas(i), temstart)  !calculate temperature
         nh(i) = rhoH(i) * dom
      end if
      enddo

!     Correct temperature for gamma from H2

      if (ispecies.gt.1) then
         do i = is+1, ie+1
         if (itmask(i)) then
            nH2 = 0.5d0*(H2I(i,j,k) + H2II(i,j,k))
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            if (nH2/nother.gt.1.0d-3) then
               x = 6100.d0/tgas(i) ! not quite self-consistent
               if (x.gt.10.d0) then
                  gamma2 = 0.5d0*5.d0
               else
                  gamma2 = 0.5d0*(5.d0 + 2.d0*x**2 * 
     &                 exp(x)/(exp(x)-1)**2)
               endif
            else
               gamma2 = 2.5d0
            endif
            gamma2 = 1.d0 + (nH2 + nother)/
     &                      (nH2*gamma2 + nother/(gamma-1.d0))
            tgas(i) = tgas(i) * (gamma2 - 1.d0)/(gamma - 1.d0)
            end if
          enddo
      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo
      return
      end
