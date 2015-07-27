subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  !use cooling_mod, only : do_radtrans, do_cool,chemistry
  use cooling_mod
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info,isink
  integer,dimension(1:nvector) :: ind_grid

  if(.not. cooling) return

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term by vector sweeps
  if (do_cool) then
     ncache=active(ilevel)%ngrid
     !$omp parallel do private(igrid,ngrid,i,ind_grid) schedule(dynamic,1)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call coolfine1(ind_grid,ngrid,ilevel)
     end do
  end if

  if (do_radtrans) then
    call radiative_cooling_fine(ilevel)
  endif
111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  use cooling_mod
  use krome_main !mandatory
  use krome_user !array sizes and utils
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer       :: i,ind,iskip,idim,nleaf
  integer, parameter :: neul=5
  real(dp)      :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)  :: dtcool
  integer,dimension(1:nvector),      save :: ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector), save :: nH,T2,delta_T2,ekk,emag
  real(kind=8), save :: time_old=-1.
  integer, save :: nprint=20

  !KROME: additional variables requested by KROME
  real*8::unoneq(krome_nmols), Tgas
  real*8::mu_noneq,mu_noneq_old,iscale_d,T2old,t2gas
  !$omp threadprivate(nH,T2,delta_T2,ekk,emag,ind_cell,ind_leaf)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),neul)
     end do

     do i=1,nleaf
        ekk(i)=0.5_8*uold(ind_leaf(i),1+1)**2/nH(i)
     end do
     do idim=2,3
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5_8*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        emag(i)=0.125_8*(uold(ind_leaf(i),1+neul)+uold(ind_leaf(i),1+nvar))**2
     end do
     do idim=2,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125_8*(uold(ind_leaf(i),idim+neul) &
                +uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
     do i=1,nleaf
        T2(i)=(uold(ind_leaf(i),ichem)-1.0)*(T2(i)-ekk(i)-emag(i)) ! Use effective adiabatic index
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute cooling time step in seconds
     dtcool = dtnew(ilevel)*scale_t

     !************************
     !KROME: inside this cooling+chemistry
     if(do_cool.or.chemistry) then

        !scale_d inverse
        iscale_d = 1d0/scale_d

        ! Compute net cooling at constant nH (original cooling)
        !call simple_cooling(nH,T2,dtcool,scale_t,delta_T2,nleaf)
        do i=1,nleaf

           ! KROME: from 2dim array of RAMSES to 1dim array for KROME
#KROME_update_unoneq

           ! KROME: from code units to 1/cm3 for KROME
#KROME_scale_unoneq

           !get the mean molecular weight
           mu_noneq_old = krome_get_mu(unoneq(:))

           !convert to K
           Tgas  = T2(i) * mu_noneq_old
           T2old = T2(i)

           ! Normalise to Av = 1 for n ~ 1e3, and let it scale like 2/3 power.
           ! This is roughly correct according to Glover et al (astro-ph:1403.3530)
           call krome_set_user_Av( (Av_rho * nH(i) / mu_noneq_old)**0.66667_8 )

           if(do_cool.and.chemistry) then
              !KROME: do chemistry+cooling
              if (any(unoneq < 0.0_dp)) then
                 write(*,*) 'Negative densities are not allowed'
                 write(*,*) 'N: ', unoneq
                 stop
              end if
              call krome(unoneq(:), Tgas, dtcool)
           elseif(.not.chemistry.and.do_cool) then
              !KROME: cooling only
              call krome_thermo(unoneq(:), Tgas, dtcool)
           elseif(.not.do_cool.and.chemistry) then
              print *,"ERROR (KROME): you cannot do chemistry without cooling"
           else
              continue
           end if

           !KROME: store the adiabatic index as the first element
           ! of the chem array (index=ichem)
           uold(ind_leaf(i),ichem) = krome_get_gamma(unoneq(:),Tgas)

           !KROME: compute mu with the chemistry updated
           mu_noneq = krome_get_mu(unoneq(:))

           !KROME: from KROME 1/cm3 to code units of RAMSES
#KROME_backscale_unoneq

           !KROME: from 1dim array of KROME to 2dim array of RAMSES
           ! indexes are shifted by 1 because of adiabatic index 
           ! position (first species index = ichem+1)
#KROME_backupdate_unoneq

           !KROME: compute t2emperature difference
           t2gas    = Tgas / mu_noneq
           delta_T2(i) = t2gas - t2old
           T2(i) = t2gas
           if (c_verbose > 0) then
           !$omp critical
           if (t .ne. time_old) then
             time_old = t
             print '(a,F10.2,F8.0,3e11.3)', 'Time [kyr], temperature [K], density [N_H cm^-3], cooling rate [erg / (s cm^3)] :', &
               t*scale_t / (krome_seconds_per_year*1e3), t2old * mu_noneq_old, nH(i),  &
               delta_T2(i)*nH(i)/scale_T2/(uold(ind_leaf(i),ichem)-1.0) * scale_d * scale_V**2 / dtcool, &
               delta_T2(i)/nH(i)/scale_T2/(uold(ind_leaf(i),ichem)-1.0) * scale_d * scale_V**2 / dtcool
             print '(a,F12.4, e12.3, 3F10.2, 12e12.3)', &
               'Time, density, T, mu, gamma, H, e-, H+, H-, H2, C, C+, O, O+, CO, OH, H2O :', &
               t * scale_t / (krome_seconds_per_year*1e3), nH(i),  &
               Tgas, mu_noneq, uold(ind_leaf(i),ichem), &
               uold(ind_leaf(i),ichem+krome_idx_H)   / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_e)   / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_Hj)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_Hk)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_H2)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_C)   / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_Cj)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_O)   / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_Oj)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_CO)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_OH)  / uold(ind_leaf(i),1), &
               uold(ind_leaf(i),ichem+krome_idx_H2O) / uold(ind_leaf(i),1)
           end if
           !$omp end critical
           end if
        end do
     end if

     ! Compute rho, new internal energy, and update total fluid energy
     do i=1,nleaf
       nH(i) = nH(i)/scale_nH
       T2(i) = T2(i)*nH(i)/scale_T2/(uold(ind_leaf(i),ichem)-1.0)
       uold(ind_leaf(i),neul) = T2(i)+ekk(i)+emag(i)
       if (T2(i) < 0.) then
         !$omp critical
         if (nprint>0) then
           nprint=nprint-1
           print '(a,i6,i4,i10,1p,6e11.3)', 'Negative temperature detected, mycpu, i, idx, Enew, rho :', &
               myid, i, ind_leaf(i), uold(ind_leaf(i),neul)-(ekk(i)+emag(i)),nH(i)
         endif
         !$omp end critical
       endif
     end do
    
  end do
  ! End loop over cells

end subroutine coolfine1

