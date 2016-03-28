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
        if (do_oct_chemistry) then
           call coolfine1_oct(ind_grid,ngrid,ilevel)
        else
           call coolfine1(ind_grid,ngrid,ilevel)
        end if
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
  real*8::mu_noneq,mu_noneq_old,iscale_d,t2old,t2gas
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
           t2old = T2(i)

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

           !KROME: compute temperature difference
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
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1_oct(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_mod
  use krome_main !mandatory
  use krome_user !array sizes and utils
  implicit none
  integer, intent(in) :: ilevel, ngrid
  integer, dimension(1:nvector), intent(in) :: ind_grid
  !-------------------------------------------------------------------
  integer :: i,ind,igrid,iskip,nvec
  integer, parameter :: neul=5
  real(kind=8)  :: dtlevel, ekin_mag, frac, Tgas
  real(kind=8), dimension(1:krome_natoms),     save :: ref
  real(kind=8), dimension(ichem_param+1:nvar), save :: x, fracV
  integer,dimension(1:nvector),                save :: nleaf
  integer,dimension(twotondim,1:nvector),      save :: ind_leaf
  ! uin and uout contains [nvector entries]x[rho, Pressure, rhoX]
  real(kind=8), dimension(1:nvector,nvar-ichem_param+2), save :: uin=1., uout=1.
  !$omp threadprivate(ref,x,fracV,nleaf,ind_leaf,uin,uout)

  ! units only needed to be computed in first call
  logical,  save :: first_call=.true.
  real(dp), save :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !$omp threadprivate(first_call,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v)

  ! Conversion factor from user units to cgs units
  if (first_call) then
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    first_call = .false.
  endif
  
  ! Compute cooling time step in seconds
  dtlevel = dtnew(ilevel) * scale_t

  ! Count leaf cells in octs
  ! ------------------------------------------------------------------------
  nleaf=0
  do ind=1,twotondim
     do igrid=1,ngrid
        iskip=ncoarse+(ind-1)*ngridmax
        if (son(ind_grid(igrid)+iskip)==0) then
           nleaf(igrid) = nleaf(igrid) + 1
           ind_leaf(nleaf(igrid),igrid) = ind_grid(igrid)+iskip
        endif
     end do
  end do

  ! Construct average states. Average over pressure instead of total energy
  ! If based on (internal) energy, the adiabatic index has to be determined implicitly
  ! ------------------------------------------------------------------------
  nvec=0
  do igrid=1,ngrid
     if (nleaf(igrid)>0) then
        nvec = nvec + 1
        uin(nvec,:)=0.0_8
        do ind=1,nleaf(igrid)
           i = ind_leaf(ind,igrid)
           uin(nvec,1) = uin(nvec,1) + uold(i,1)
           ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
                      0.125_8*( &
                        (uold(i,1+neul) + uold(i,1+nvar))**2 + &
                        (uold(i,2+neul) + uold(i,2+nvar))**2 + &
                        (uold(i,3+neul) + uold(i,3+nvar))**2 )
           uin(nvec,2)  = uin(nvec,2) + (uold(i,neul) - ekin_mag) * (uold(i,ichem)-1.0_8)
           uin(nvec,3:) = uin(nvec,3:) + uold(i,ichem+1:nvar)
        end do
        uin(nvec,:) = uin(nvec,:) / nleaf(igrid)
     end if
  end do

  ! Compute chemical evolution of oct-averaged states over time dtcool
  ! State vector contains [rho, pressure, rhoX]
  ! ------------------------------------------------------------------------
  call evolve_chemistry_and_Tgas(uin,nvec,dtlevel,uout)

  ! Scatter results back to cells
  ! ------------------------------------------------------------------------
  nvec=0
  do igrid=1,ngrid
     ! Only one leaf cell in oct. Copy over results
     if (nleaf(igrid)==1) then
        nvec  = nvec + 1
        x = uout(nvec,3:)

        Tgas = uout(nvec,2) / uout(nvec,1) * scale_T2 * krome_get_mu_x(x)   ! Compute temperature from pressure
        if (Tgas < 0.) then
           !$omp critical
           print '(a,i6,i10,i3,1p,40e11.3)', 'Negative temperature detected, mycpu, i, idx, Tgas, rho :', &
             myid, i, ind_leaf(ind,igrid), Tgas, uout(nvec,1), x
           stop
           !$omp end critical
        endif
  
        i = ind_leaf(1,igrid)
        uold(i,ichem) = krome_get_gamma_x(x,Tgas)                           ! New adiabatic index 
        uold(i,ichem+1:nvar) = x                                            ! Update abundances
        ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
                      0.125_8*( &
                           (uold(i,1+neul) + uold(i,1+nvar))**2 + &
                           (uold(i,2+neul) + uold(i,2+nvar))**2 + &
                           (uold(i,3+neul) + uold(i,3+nvar))**2 )
        uold(i,neul) = uout(2,nvec) / (uold(i,ichem) - 1.0_8) + ekin_mag    ! Total energy
     endif
     ! More than one leaf cell. Rescale fractionally, and make sure we conserve metallicities
     if (nleaf(igrid)>1) then
        nvec  = nvec + 1

        frac = uout(nvec,2) / uin(nvec,2)
        fracV = uout(nvec,3:) / uin(nvec,3:)                                ! fractional change in species
        do ind=1,nleaf(igrid)
           i = ind_leaf(ind,igrid)

           ! Update chemical species according to fractional change
           ! Make sure interpolated chemical densities conserve metallicity
           ! ---------------------------------------------------------------
           x = uold(i,ichem+1:nvar)                                         ! extract current abundances
           ref = krome_conserveLinGetRef_x(x)                               ! get reference metallicity
           x = x * fracV                                                    ! rescale abundances 
           call krome_conserveLin_x(x,ref)                                  ! make sure metallicity is conserved

           uold(i,ichem+1:nvar) = x                                         ! store updated abundance

           ! Update pressure according to fractional change
           ! Compute Tgas and find adiabatic index. Find new total energy
           ! ---------------------------------------------------------------
           ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
                      0.125_8*( &
                        (uold(i,1+neul) + uold(i,1+nvar))**2 + &
                        (uold(i,2+neul) + uold(i,2+nvar))**2 + &
                        (uold(i,3+neul) + uold(i,3+nvar))**2 )
           uold(i,neul) = (uold(i,ichem) - 1.0_8) * (uold(i,neul) - ekin_mag)*frac ! Update pressure
           
           Tgas = uold(i,neul) / uold(i,1) * scale_T2 * krome_get_mu_x(x)   ! Compute temperature from pressure

           if (Tgas < 0.) then
              !$omp critical
              print '(a,i6,i10,i3,1p,40e11.3)', 'Negative temperature detected, mycpu, i, idx, Tgas, rho :', &
                myid, i, ind_leaf(ind,igrid), Tgas, uold(i,1), uold(i,ichem+1:nvar)
              stop
              !$omp end critical
           endif

           uold(i,ichem) = krome_get_gamma_x(x,Tgas)                        ! New adiabatic index 
           uold(i,neul) =  uold(i,neul) / (uold(i,ichem) - 1.0_8) + ekin_mag! New total energy
        end do
     end if
  end do

end subroutine coolfine1_oct
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine evolve_chemistry_and_Tgas(uin,nvec,dt,uout)
  use amr_commons
  use hydro_commons
  use cooling_module
  use cooling_mod
  use krome_main !mandatory
  use krome_user !array sizes and utils
  implicit none
  integer, intent(in) :: nvec
  real(kind=8), dimension(nvector,nvar-ichem_param+2), intent(in) :: uin
  real(kind=8), dimension(nvector,nvar-ichem_param+2), intent(out):: uout
  real(kind=8), intent(in) :: dt
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer       :: i
  integer, parameter :: neul=5
  real(kind=8), dimension(nvector),     save :: nH,T2
  real(kind=8), dimension(krome_nmols), save :: unoneq
  real(kind=8) :: Tgas, mu_noneq,mu_noneq_old,iscale_d,t2gas
  !$omp threadprivate(nH,T2,unoneq)

  ! units only needed to be computed in first call
  logical,  save :: first_call=.true.
  real(dp), save :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !$omp threadprivate(first_call,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v)

  ! Conversion factor from user units to cgs units
  if (first_call) then
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    first_call = .false.
  endif
  
  ! Compute rho
  nH = uin(:,1)

  ! Compute T2=T/mu in Kelvin
  T2 = uin(:,2) / nH * scale_T2

  ! Compute nH in H/cc
  nH = nH * scale_nH

  !scale_d inverse
  iscale_d = 1.0_8/scale_d

  ! Compute net cooling at constant nH (original cooling)
  do i=1,nvec

    ! KROME: from 2dim array of RAMSES to 1dim array for KROME
#KROME_vecupdate_unoneq

    ! KROME: from code units to 1/cm3 for KROME
#KROME_scale_unoneq

    !get the mean molecular weight
    mu_noneq_old = krome_get_mu(unoneq(:))

    !convert to K
    Tgas  = T2(i) * mu_noneq_old

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
      call krome(unoneq(:), Tgas, dt)
    elseif(.not.chemistry.and.do_cool) then
      !KROME: cooling only
      call krome_thermo(unoneq(:), Tgas, dt)
    elseif(.not.do_cool.and.chemistry) then
      write(*,*) 'ERROR (KROME): you cannot do chemistry without cooling'
    else
      continue
    end if

    !KROME: compute mu with the chemistry updated
    mu_noneq = krome_get_mu(unoneq(:))

    !KROME: from KROME 1/cm3 to code units of RAMSES
#KROME_backscale_unoneq

    !KROME: from 1dim array of KROME to 2dim array of RAMSES
    ! indexes are shifted by 1 because of adiabatic index 
    ! position (first species index = ichem+1)
#KROME_vecbackupdate_unoneq

    !KROME: compute temperature difference
    T2(i) = Tgas / mu_noneq

    if (Tgas < 0.) then
      !$omp critical
      print '(a,i6,i3,1p,2e11.3)', 'Negative temperature detected, mycpu, i, T, rho :', myid, i, Tgas, nH(i)
      stop
      !$omp end critical
    endif

  end do
  ! End loop over cells

  ! Compute new pressure
  uout(1:nvec,2) = T2(1:nvec) * uin(1:nvec,1) / scale_T2

end subroutine evolve_chemistry_and_Tgas
