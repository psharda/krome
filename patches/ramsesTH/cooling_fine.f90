subroutine cooling_fine(ilevel)
  use amr_parameters, only : ifout
  use amr_commons
  use hydro_commons
  use cooling_module
  !use cooling_mod, only : do_radtrans, do_cool,chemistry
  use cooling_mod
  use ray_utils_m, only : file_location
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
  character(len=80) :: filename

  integer, parameter :: u_dbg = 204   ! File unit for debug dump
  logical :: dump_temperature
  character(len=255) ::dbg_fname     ! File name for debug dump

  if(.not. cooling) return

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  if (do_radtrans) then
    if(c_verbose > 0 .and. myid==1) print*, "Solving radiative transfer"
    call timer("hydro - cooling - lampray", "start")
    call make_virtual_fine_dp(uold(1,1),ilevel,3+nvar)
    call radiative_cooling_fine(ilevel)
    call timer('hydro - cooling','start)')
  endif

  ! Open file for debug dump of temperature
  dump_temperature = (c_verbose > 1)
  if(dump_temperature) then
     write(dbg_fname, "(a,i5.5,a)") "temperature_",myid,".dat"
     open(u_dbg, file=trim(dbg_fname))
  end if
  
  ! Operator splitting step for cooling source term by vector sweeps
  if (do_cool) then
    if(c_verbose > 0 .and. myid==1) print*, "Solving chemistry with KROME"
    call timer("hydro - cooling - KROME", 'start') 
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
    call timer('hydro - cooling','start')
 end if

   ! Close file for debug dump of temperature
  if(dump_temperature) then
     close(u_dbg)
  end if


111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_parameters, only : icoarse_min, jcoarse_min, kcoarse_min
  use amr_commons
  use hydro_commons
  use cooling_module
  use cooling_mod
  use krome_main !mandatory
  use krome_user !array sizes and utils
  use lampray, only : radiation_field
  use lampray_parameters, only : nBin, nBinTot, nBinSS
#ifdef SELFSHIELDING
  use lampray_parameters, only : ibin_H2, ibin_CO
  use richings_dissociation_rates, only : gamma_H2_thin, gamma_CO_thin
#endif

  use stars_m,          only : xc
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer       :: i,ind,iskip,idim,nleaf
  integer, parameter :: neul=5
  real(dp)      :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)  :: dtcool
  integer,dimension(1:nvector),      save :: ind_cell,ind_leaf,ind_grid_leaf
  real(kind=8),dimension(1:nvector), save :: nH,T2,delta_T2,ekk,emag
  real(kind=8),dimension(1:nvector,3), save :: xleaf
  real(kind=8), save :: time_old=-1.
  integer, save :: nprint=20
  real*8::phbin(nBin)
  real*8::phbin_ss(nBinSS)
  real(kind=8) :: T_tmp
  real(kind=8) :: c(3)

  !KROME: additional variables requested by KROME
  real*8::unoneq(krome_nmols), Tgas
  real*8::mu_noneq,mu_noneq_old,iscale_d,t2old,t2gas

  integer, parameter :: u_dbg = 204   ! File unit for debug dump
  logical :: dump_temperature
  
  !$omp threadprivate(nH,T2,delta_T2,ekk,emag,ind_cell,ind_leaf,ind_grid_leaf,xleaf)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  dump_temperature = (c_verbose > 1)

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
        ind_grid_leaf(nleaf)=ind_grid(i)
        ind_leaf(nleaf)=ind_cell(i)
      end if
    end do

    if(c_verbose > 1) then 
      ! Store leaf position for debug output
      do i=1,nleaf
        xleaf(i,:) = xg(ind_grid_leaf(i),:) + xc(ind,:,ilevel)
      end do
    end if

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

        if(do_radtrans) then
          ! Set intensity from RT
          phbin = radiation_field(ind_leaf(i),1:nBin)
          if(any(phbin.lt.0.0_dp)) then
             write(*,*) 'Warning: Negative intensity found in cell ', ind_leaf(i)
             phbin = merge(phbin, 0.0_dp, phbin.ge.0.0_dp)
             !call clean_stop
          end if
          call krome_set_photoBinJ(phbin)
          phbin_ss = radiation_field(ind_leaf(i),nBin+1:nBinTot)
#ifdef SELFSHIELDING
          call krome_set_user_gamma_H2(phbin_ss(ibin_H2)*gamma_H2_thin)
          call krome_set_user_gamma_CO(phbin_ss(ibin_CO)*gamma_CO_thin)
#endif
        else
          ! Normalise to Av = 1 for n ~ 1e3, and let it scale like 2/3 power.
          ! This is roughly correct according to Glover et al (astro-ph:1403.3530)
          call krome_set_user_Av( (Av_rho * nH(i) / mu_noneq_old)**0.66667_8 )
        end if

        if(do_cool.and.chemistry) then
          !KROME: do chemistry+cooling
          if (any(unoneq < 0.0_dp)) then
            write(*,*) 'Negative densities are not allowed'
            write(*,*) 'N: ', unoneq
            stop
         end if

         c(1)  = xg(ind_grid_leaf(i),1) + xc(ind,1,ilevel) - icoarse_min
#if NDIM>1
         c(2)  = xg(ind_grid_leaf(i),2) + xc(ind,2,ilevel) - jcoarse_min
#endif
#if NDIM>2
         c(3)  = xg(ind_grid_leaf(i),3) + xc(ind,3,ilevel) - kcoarse_min
#endif
         if(dump_temperature) then
            !$omp critical
            write(u_dbg,*) c, Tgas
            !$omp end critical
         endif
         
          call krome(unoneq(:), Tgas, dtcool)
        elseif(.not.chemistry.and.do_cool) then
          !KROME: cooling only
          call krome_thermo(unoneq(:), Tgas, dtcool)
        elseif(.not.do_cool.and.chemistry) then
          print *,"ERROR (KROME): you cannot do chemistry without cooling"
          stop
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
                !!$             print '(a,F12.4, e12.3, 3F10.2, 12e12.3)', &
                !!$               'Time, density, T, mu, gamma, H, e-, H+, H-, H2, C, C+, O, O+, CO, OH, H2O :', &
                !!$               t * scale_t / (krome_seconds_per_year*1e3), nH(i),  &
                !!$               Tgas, mu_noneq, uold(ind_leaf(i),ichem), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_H)   / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_e)   / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_Hj)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_Hk)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_H2)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_C)   / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_Cj)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_O)   / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_Oj)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_CO)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_OH)  / uold(ind_leaf(i),1), &
                !!$               uold(ind_leaf(i),ichem+krome_idx_H2O) / uold(ind_leaf(i),1)
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
  integer, save :: noct_evolved=0, nrollback=0, ncell_evolved=0, nprint_threshold=1000
  integer, parameter :: neul=5
  real(kind=8)  :: dtlevel, ekin_mag, frac, Tgas
  real(kind=8), dimension(1:krome_natoms),     save :: ref
  real(kind=8), dimension(ichem_param+1:nvar), save :: x, xout, fracV
  real(kind=8), dimension(ichem_param+1:nvar,twotondim), save :: xvec
  integer,dimension(1:nvector),                save :: nleaf
  integer,dimension(twotondim,1:nvector),      save :: ind_leaf
  logical, external :: similar_states
  logical :: no_rollback
  ! uin and uout contains [nvector entries]x[rho, Pressure, rhoX]
  integer, parameter :: nstate = max(nvector,twotondim)
  real(kind=8), dimension(nstate,nvar-ichem_param+2), save :: uin=1., uout=1.
  !$omp threadprivate(ref,x,xout,xvec,fracV,nleaf,ind_leaf,uin,uout)

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
      if (similar_states(ind_leaf(:,igrid),nleaf(igrid))) then
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
      else
        call solve_chemistry_in_oct(ind_leaf(:,igrid),nleaf(igrid),dtlevel)
        !$omp atomic
        ncell_evolved = ncell_evolved + nleaf(igrid)
      endif
    endif
  enddo

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
      uold(i,neul) = uout(nvec,2) / (uold(i,ichem) - 1.0_8) + ekin_mag    ! Total energy
      !$omp atomic
      ncell_evolved = ncell_evolved + nleaf(igrid)
    endif
    ! More than one leaf cell. Rescale fractionally, and make sure we conserve metallicities
    if (nleaf(igrid)>1) then
      nvec  = nvec + 1

      if (any(uin(nvec,:).ne.uin(nvec,:))) then
        write(*,*) 'Problems with input state vector :'
        write(*,*) uin(nvec,:)
        write(*,*) uout(nvec,:)
        write(*,*) uold(ind_leaf(1,igrid),1),uold(ind_leaf(1,igrid),ichem+1:nvar)
        stop
      endif
      frac = uout(nvec,2) / uin(nvec,2)
      fracV = uout(nvec,3:) / (uin(nvec,3:) + tiny(1.0_8))                ! fractional change in species

      ! compute cell-by-cell updated abundances and check if we need serious corrections
      no_rollback=.true.
      do ind=1,nleaf(igrid)
        i = ind_leaf(ind,igrid)

        ! Update chemical species according to fractional change
        ! Make sure interpolated chemical densities conserve metallicity
        ! ---------------------------------------------------------------
        xout   = uold(i,ichem+1:nvar)                                    ! extract current abundances
        ref = krome_conserveLinGetRef_x(xout)                            ! get reference metallicity
        xout = xout * fracV                                              ! rescale abundances
        xvec(:,ind) = xout
        call krome_conserveLin_x(xvec(:,ind),ref)                        ! make sure metallicity is conserved

        ! If there is more than a 10% change in any significant species then
        ! roll back and compute chemistry cell-by-cell in oct
        if (any(abs(xout - xvec(:,ind)) / (xout + 1d-10*uold(i,1)) > 0.1)) then
          call solve_chemistry_in_oct(ind_leaf(:,igrid), nleaf(igrid),dtlevel)
          no_rollback=.false.
          !$omp atomic
          nrollback = nrollback + nleaf(igrid)
          exit                                                                  ! jump out of loop over ind
        endif
      enddo

      ! scatter results to cells
      if (no_rollback) then
        !$omp atomic
        noct_evolved = noct_evolved + nleaf(igrid)
        do ind=1,nleaf(igrid)
          i = ind_leaf(ind,igrid)
          ! Update pressure according to fractional change
          ! Compute Tgas and find adiabatic index. Find new total energy
          ! ---------------------------------------------------------------
          ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
              0.125_8*( &
              (uold(i,1+neul) + uold(i,1+nvar))**2 + &
              (uold(i,2+neul) + uold(i,2+nvar))**2 + &
              (uold(i,3+neul) + uold(i,3+nvar))**2 )
          uold(i,neul) = (uold(i,ichem) - 1.0_8) * (uold(i,neul) - ekin_mag)*frac ! Update pressure

          Tgas = uold(i,neul) / uold(i,1) * scale_T2 * krome_get_mu_x(xvec(:,ind))   ! Compute temperature from pressure

          if (Tgas < 0.) then
            !$omp critical
            print '(a,i6,i10,i3,1p,40e11.3)', 'Negative temperature detected, mycpu, i, idx, Tgas, rho :', &
                myid, i, ind_leaf(ind,igrid), Tgas, uold(i,1), uold(i,ichem+1:nvar)
            stop
            !$omp end critical
          endif

          uold(i,ichem) = krome_get_gamma_x(xvec(:,ind),Tgas)              ! New adiabatic index
          uold(i,ichem+1:nvar) = xvec(:,ind)                               ! Update abundances
          uold(i,neul) =  uold(i,neul) / (uold(i,ichem) - 1.0_8) + ekin_mag! New total energy
        enddo
      endif
    endif
  enddo

  !$omp critical
  if ((nrollback + noct_evolved + ncell_evolved) > nprint_threshold) then
    write(*,'(a,i4.4,3(a,i9),a,3(1x,i3),a)') 'C=',myid,'coolfine1_oct: noct_evolved=',noct_evolved, &
        ' ncell_evolved=', ncell_evolved, ' nrollback=', nrollback, &
        ' Fraction=', &
        (100*noct_evolved)  / (nrollback + noct_evolved + ncell_evolved), &
        (100*ncell_evolved) / (nrollback + noct_evolved + ncell_evolved), &
        (100*nrollback)     / (nrollback + noct_evolved + ncell_evolved), &
        '%'
    nprint_threshold = nprint_threshold * 1.5
  endif
  !$omp end critical
end subroutine coolfine1_oct
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function similar_states(ind_leaf,nleaf)
  use amr_commons
  use hydro_commons
  implicit none
  integer, dimension(twotondim), intent(in) :: ind_leaf
  integer,                       intent(in) :: nleaf
  logical                                   :: similar_states
  !
  integer, parameter :: neul=5
  integer :: i, isp, ind, idx
  real(kind=8) :: x, av, mn, mx, rms, iav, eps, ekin_mag
  !
  real(kind=8), parameter :: tol1 = 0.05 !0.4  ! Allow 40% relative difference between average and min or max
  real(kind=8), parameter :: tol2 = 2./3.*tol1 ! Allow 27% relative difference for rms (tol1+tol2 consistent for normal distributed data)
  !
  similar_states = .true.
  do isp=ichem_param-2,nvar
    av=0.0_8
    rms = 0.0_8
    mn=huge(1.0_8)
    mx=0.0_8
    idx = isp
    if (isp==ichem_param-2) idx=1 ! density
    if (isp==ichem_param-1) idx=5 ! pressure
    !
    ! calculate average, rms, minimum, and maximum for species
    !
    do ind=1,nleaf
      i = ind_leaf(ind)
      x = uold(i,idx)
      if (idx==5) then
        ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
            0.125_8*( &
            (uold(i,1+neul) + uold(i,1+nvar))**2 + &
            (uold(i,2+neul) + uold(i,2+nvar))**2 + &
            (uold(i,3+neul) + uold(i,3+nvar))**2 )
        x  = (x - ekin_mag) * (uold(i,ichem)-1.0_8)
      endif
      av = av + x
      rms = rms + x*x
      mn = min(mn,x)
      mx = max(mx,x)
    enddo
    av = av / nleaf
    rms = sqrt(abs(rms / nleaf - av*av)) ! avoid round-off making it negative

    mn = mn + tiny(1.0_8)        ! guard for 0 / 0
    av = av + tiny(1.0_8)        ! guard for 0 / 0
    mx = mx + tiny(1.0_8)        ! guard for 0 / 0

    if (idx==1) eps = mn * 1d-10 ! The tolerance for disregarding species is at 10^-10 of the minimum density in the cells

    ! check if state vector is too disimilar for density, pressure, or chemical mass fraction
    iav = 1.0_8 / (av + eps)
    if ((av-mn)*iav > tol1 .or. (mx-av)*iav > tol1 .or. rms*iav > tol2) then
      similar_states = .false.
      exit
    endif
  enddo
end function similar_states
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine solve_chemistry_in_oct(ind_leaf,nleaf,dtlevel)
  use amr_commons
  use hydro_commons
  use cooling_mod
  use krome_main !mandatory
  use krome_user !array sizes and utils
  implicit none
  integer, dimension(twotondim), intent(in) :: ind_leaf
  integer,                       intent(in) :: nleaf
  real(kind=8),                  intent(in) :: dtlevel
  !
  integer, parameter :: neul=5
  integer            :: i, ind
  real(kind=8)       :: ekin_mag, Tgas
  real(kind=8), dimension(ichem_param+1:nvar), save :: x
  ! uin and uout contains [2^ndim entries]x[rho, Pressure, rhoX]
  integer, parameter :: nstate = max(nvector,twotondim)
  real(kind=8), dimension(nstate,nvar-ichem_param+2), save :: uin, uout
  !$omp threadprivate(x,uin,uout)
  !
  do ind=1,nleaf
    i = ind_leaf(ind)
    uin(ind,1) = uold(i,1)
    ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
        0.125_8*( &
        (uold(i,1+neul) + uold(i,1+nvar))**2 + &
        (uold(i,2+neul) + uold(i,2+nvar))**2 + &
        (uold(i,3+neul) + uold(i,3+nvar))**2 )
    uin(ind,2)  = (uold(i,neul) - ekin_mag) * (uold(i,ichem)-1.0_8)
    uin(ind,3:) = uold(i,ichem+1:nvar)
  end do

  ! Compute chemical evolution of cell states over time dtlevel
  ! State vector contains [rho, pressure, rhoX]
  ! ------------------------------------------------------------------------
  call evolve_chemistry_and_Tgas(uin,nleaf,dtlevel,uout)

  do ind=1,nleaf
    x = uout(ind,3:)
    if (abs(sum(x) - uout(ind,1)) > 1d-5*(uout(ind,1)+tiny(1.0_8)) .or. uout(ind,1) < 1e-100_8) then
      write(*,*) 'Species do not sum up to total density or density is zero :', sum(x), uin(ind,1), uout(ind,1)
      write(*,*) 'Input state vector  :', uin(ind,:)
      write(*,*) 'Output state vector :', uout(ind,:)
      stop
    endif
    Tgas = uout(ind,2) / uout(ind,1) * krome_get_mu_x(x)   ! Compute temperature from pressure
    if (Tgas < 0.) then
      !$omp critical
      print '(a,i6,i10,i3,1p,40e11.3)', 'Negative temperature detected, mycpu, i, idx, Tgas (numerical units), rho :', &
          myid, i, ind_leaf(ind), Tgas, uout(ind,1), x
      stop
      !$omp end critical
    endif
    i = ind_leaf(ind)
    uold(i,ichem) = krome_get_gamma_x(x,Tgas)                           ! New adiabatic index
    uold(i,ichem+1:nvar) = x                                            ! Update abundances
    ekin_mag = 0.5_8 * (uold(i,2)**2 + uold(i,3)**2 + uold(i,4)**2) / uold(i,1) + &
        0.125_8*( &
        (uold(i,1+neul) + uold(i,1+nvar))**2 + &
        (uold(i,2+neul) + uold(i,2+nvar))**2 + &
        (uold(i,3+neul) + uold(i,3+nvar))**2 )
    uold(i,neul) = uout(ind,2) / (uold(i,ichem) - 1.0_8) + ekin_mag    ! Total energy
    uold(i,ichem+1:nvar) = x                                            ! store updated abundance
  enddo
end subroutine solve_chemistry_in_oct
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
  integer, parameter :: nstate = max(nvector,twotondim)
  real(kind=8), dimension(nstate,nvar-ichem_param+2), intent(in) :: uin
  real(kind=8), dimension(nstate,nvar-ichem_param+2), intent(out):: uout
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
  nH(1:nvec) = uin(1:nvec,1)

  ! Compute T2=T/mu in Kelvin
  T2(1:nvec) = uin(1:nvec,2) / nH(1:nvec) * scale_T2

  ! Compute nH in H/cc
  nH(1:nvec) = nH(1:nvec) * scale_nH

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

  ! Set total density
  uout(1:nvec,1) = uin(1:nvec,1)

  ! Compute new pressure
  uout(1:nvec,2) = T2(1:nvec) * uin(1:nvec,1) / scale_T2

end subroutine evolve_chemistry_and_Tgas
