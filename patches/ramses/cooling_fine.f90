subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Compute sink accretion rates
  if(sink)call compute_accretion_rate(0)

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do
  
  if(cooling.and.ilevel==levelmin.and.cosmo.and..not.chemistry)then
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
  end if

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use hydro_parameters, ONLY: nvar
  use krome_main
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz,indchem
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant,rhogas,mu_noneqold,mu_noneq
  real(dp)::dthydro
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::nH,T2,delta_T2,ekk
  real(kind=8),dimension(1:nvector),save::t2gas,t2gasold,tgas,tgasold
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(kind=8),dimension(1:11),save::unoneq
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  endif

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
     
     ! Compute metallicity in solar units
     if(metal)then
        do i=1,nleaf
           Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
        end do
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

     ! Compute pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=exp(-nH(i)/0.01)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     ! Compute net cooling at constant density
     if(chemistry)then
        do i=1,nleaf
           ! Compute "thermal" temperature by substracting polytrope
           T2(i) = max(T2(i)-T2min(i),T2_min_fix)
	   ! Species mass density in code units
           unoneq(1) = uold(ind_leaf(i),ndim+4) ! HI                    
           unoneq(2) = uold(ind_leaf(i),ndim+5) ! E  
           unoneq(3) = uold(ind_leaf(i),ndim+6) ! HII  
           unoneq(4) = uold(ind_leaf(i),ndim+7) ! HeI  
           unoneq(5) = uold(ind_leaf(i),ndim+8) ! HeII 
           unoneq(6) = uold(ind_leaf(i),ndim+9) ! HeIII
           unoneq(7) = uold(ind_leaf(i),ndim+10)! HM
           unoneq(8) = uold(ind_leaf(i),ndim+11)! H2I  
           unoneq(9) = uold(ind_leaf(i),ndim+12)! H2II 
           unoneq(10)= 0.0d0                    ! CR 
           unoneq(11)= 0.0d0                    ! g            
           rhogas     = (nH(i)/scale_nH)        ! code units
           ! Compute the non-eq. mean molecular weight
	   call get_mu(unoneq(:),rhogas,mu_noneqold)
           t2gasold(i)= T2(i)               
           tgasold(i) = T2(i)*mu_noneqold       ! K
           t2gas(i)   = t2gasold(i)             ! K
           tgas(i)    = tgasold(i)
           dthydro    = dtcool                  ! s
           ! Species mass density (in code units) -> number denstiy for krome 1/cm3
           unoneq(1) = unoneq(1)*scale_d/mH        ! HI                    
           unoneq(2) = unoneq(2)*scale_d/me        ! E  
           unoneq(3) = unoneq(3)*scale_d/mH        ! HII  
           unoneq(4) = unoneq(4)*scale_d/(4.d0*mH) ! HeI  
           unoneq(5) = unoneq(5)*scale_d/(4.d0*mH) ! HeII 
           unoneq(6) = unoneq(6)*scale_d/(4.d0*mH) ! HeIII
           unoneq(7) = unoneq(7)*scale_d/mH        ! HM
           unoneq(8) = unoneq(8)*scale_d/(2.d0*mH) ! H2I  
           unoneq(9) = unoneq(9)*scale_d/(4.d0*mH) ! H2II 
           unoneq(10)= 0.0d0                       ! CR 
           unoneq(11)= 0.0d0                       ! g 
           ! Call non-equilibrium cooling package KROME
           ! with c.g.s. parameters 
           call krome(unoneq(:),tgas(i),dthydro)
           ! Species number denstiy (in c.g.s.) -> mass density (code units) 
           unoneq(1) = unoneq(1)*mH/scale_d        ! HI                    
           unoneq(2) = unoneq(2)*me/scale_d        ! E  
           unoneq(3) = unoneq(3)*mH/scale_d        ! HII  
           unoneq(4) = unoneq(4)*(4.d0*mH)/scale_d ! HeI  
           unoneq(5) = unoneq(5)*(4.d0*mH)/scale_d ! HeII 
           unoneq(6) = unoneq(6)*(4.d0*mH)/scale_d ! HeIII
           unoneq(7) = unoneq(7)*mH/scale_d        ! HM
           unoneq(8) = unoneq(8)*(2.d0*mH)/scale_d ! H2I  
           unoneq(9) = unoneq(9)*(4.d0*mH)/scale_d ! H2II 
           unoneq(10)= 0.0d0                     ! CR 
           unoneq(11)= 0.0d0                     ! g 
           ! print chemical evolution
           ! Update the species (mass density) array in code units
           do indchem=1,9
              uold(ind_leaf(i),ndim+3+indchem) = unoneq(indchem)
           end do
           ! Save gas temperature in K for the output
           uold(ind_leaf(i),ndim+3) = tgas(i)
           ! Compute new mu
	   call get_mu(unoneq(:),rhogas,mu_noneq) 
           ! Compute t2emperature difference
           t2gas(i)    = tgas(i)/mu_noneq
           delta_T2(i) = t2gas(i) - t2gasold(i)
        enddo
     else if(cooling) then
        ! Compute "thermal" temperature by substracting polytrope
        do i=1,nleaf
           T2(i) = max(T2(i)-T2min(i),T2_min_fix)
        end do
        call solve_cooling(nH,T2,Zsolar,boost,dtcool,delta_T2,nleaf)
     end if

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Compute net energy sink
     if(chemistry.or.cooling)then
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch=uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1)
              if(cooling_switch>1d-3)then
                 delta_T2(i)=0
              endif
           end do
        endif
     end if
     ! Compute minimal total energy from polytrope
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0) + ekk(i)
     end do
     
     ! Update total fluid energy
     do i=1,nleaf
        T2(i) = uold(ind_leaf(i),ndim+2)
     end do
     if(cooling.or.chemistry)then
        do i=1,nleaf
           T2(i) = T2(i)+delta_T2(i)
        end do
     endif
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i)
        end do
     else
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = max(T2(i),T2min(i))
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=20d0*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=uold(ind_leaf(i),idelay)*damp_factor
        end do
     endif
     
  end do
  ! End loop over cells

end subroutine coolfine1
!**************************************************
! Function to compute the non-eq. molecular weight
!**************************************************
subroutine get_mu(unoneq,rhogas,molweight)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp),dimension(1:11)::unoneq
  real(dp)::rhogas,molweight
  real(dp)::xHI,xE,xHII,xHeI,xHeII,xHeIII,xHM,xH2I,xH2II
  real(dp)::zHII,zHeII,zHeIII,zHM,zH2II
  real(dp)::fiHII,fiHeII,fiHeIII,fiHM,fiH2II
  real(dp)::mu_i,mu_e
  ! Species mass fractions, charges
  ! and ionization fraction to compute mu
  xHI     = unoneq(1)/rhogas ; 
  xE      = unoneq(2)/rhogas ; 
  xHII    = unoneq(3)/rhogas ; zHII   = 1.d0 ; fiHII   = xHII/(xHI+xHII+xHM+xH2I+xH2II)
  xHeI    = unoneq(4)/rhogas ; 
  xHeII   = unoneq(5)/rhogas ; zHeII  = 1.d0 ; fiHeII  = xHeII/(xHeI+xHeII+xHeIII)
  xHeIII  = unoneq(6)/rhogas ; zHeIII = 2.d0 ; fiHeIII = xHeIII/(xHeI+xHeII+xHeIII)
  xHM     = unoneq(7)/rhogas ; zHM    =-1.d0 ; fiHM    = xHM/(xHI+xHII+xHM+xH2I+xH2II) 
  xH2I    = unoneq(8)/rhogas ; 
  xH2II   = unoneq(9)/rhogas ; zH2II  = 1.d0 ; fiH2II  = xH2II/(xHI+xHII+xHM+xH2I+xH2II)
  ! Mean molecular weight
  mu_i       = xHI + xHII + xHM +xH2I/2.d0 + xH2II/2.d0 + xHeI/4.d0 + xHeII/4.d0 + xHeIII/4.d0
  mu_e       = zHII*xHII*fiHII + zHeII*(xHeII/4.d0)*fiHeII + zHeIII*(xHeIII/4.d0)*fiHeIII +zHM*(xHM)*fiHM + zH2II*(xH2II/2.d0)*fiH2II
  mu_i       = 1.d0/mu_i
  mu_e       = 1.d0/mu_e
  molweight  = 1.d0/(1.d0/mu_i + 1.d0/mu_e)

  return
end subroutine get_mu

