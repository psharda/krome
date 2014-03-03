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

  ! Compute sink accretion rates
  if(sink)call compute_accretion_rate(0)

  ! Operator splitting step for cooling source term by vector sweeps
  if (do_cool) then
     ncache=active(ilevel)%ngrid
     !!$omp parallel do private(igrid,ngrid,i,ind_grid)
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call coolfine1(ind_grid,ngrid,ilevel)
     end do
     if(ilevel==levelmin.and.cosmo)then
        if(myid==1)write(*,*)'Computing new cooling table'
        call set_table(dble(aexp))
     end if
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
  integer, save :: nprint=20

  !KROME: additional variables requested by KROME
  real*8::unoneq(krome_nmols), Tgas(nvector)
  real*8::mu_noneq,iscale_d,T2old,t2gas
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
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-emag(i))
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

           !get the mean molecular weight
           mu_noneq = krome_get_mu(unoneq(:))

           !convert to K
           Tgas(i) = T2(i) * mu_noneq
           T2old = Tgas(i)

           ! KROME: from code units to 1/cm3 for KROME
#KROME_scale_unoneq

           if(do_cool.and.chemistry) then
              !KROME: do chemistry+cooling
              call krome(unoneq(:), Tgas(i), dtcool)
           elseif(.not.chemistry.and.do_cool) then
              !KROME: cooling only
              call krome_thermo(unoneq(:), Tgas(i), dtcool)
           elseif(.not.do_cool.and.chemistry) then
              print *,"ERROR (KROME): you cannot do chemstry without cooling"
           else
              continue
           end if

           !KROME: from KROME 1/cm3 to code units of RAMSES
#KROME_backscale_unoneq

           !KROME: from 1dim array of KROME to 2dim array of RAMSES
           ! indexes are shifted by 1 because of adiabatic index 
           ! position (first species index = ichem+1)
#KROME_backupdate_unoneq

           !KROME: store the adiabatic index as the first element
           ! of the chem array (index=ichem)
           uold(ind_leaf(i),ichem) = krome_get_gamma(unoneq(:),Tgas(i))

           !Save gas temperature in K for the output
           !uold(ind_leaf(i),ndim+3) = tgas(i)

           !KROME: compute mu with the chemistry updated
           mu_noneq = krome_get_mu(unoneq(:))

           !KROME: compute t2emperature difference
           t2gas    = Tgas(i) / mu_noneq
           delta_T2(i) = t2gas - t2old
        end do
     end if

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink, 
     do i=1,nleaf
        delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
     end do
     


     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),neul)
        uold(ind_leaf(i),neul) = T2(i)+delta_T2(i)
        if (uold(ind_leaf(i),neul) < ekk(i)+emag(i)) then
           !$omp critical
           if (nprint>0) then
              nprint=nprint-1
              print '(a,i6,i4,i10,1p,6e11.3)', 'Negative energy detected, mycpu, i, idx, Eold, dT, Enew :', &
                   myid, i, ind_leaf(i), T2(i)-(ekk(i)+emag(i)), &
                   delta_T2(i), uold(ind_leaf(i),neul)-(ekk(i)+emag(i)),ekk(i),emag(i),nH(i)
           endif
           !$omp end critical
        endif
     end do

  end do
  ! End loop over cells

end subroutine coolfine1

