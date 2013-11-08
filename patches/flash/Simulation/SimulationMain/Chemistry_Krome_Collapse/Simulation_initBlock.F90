!!****if* source/Simulation/SimulationMain/Chemsitry_Krome_Collapse/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!! DESCRIPTION
!!
!! Sets up minihalo with constant density for the Krome collapse test in 1D/3D
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_initblock (blockID)



!!***used modules from FLASH.***

  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_putPointData
  
  use Eos_interface, ONLY : Eos 
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
  
  ! HERE are the arguments
  integer, intent(in) :: blockID
  
  
  
  !!***block number and (cell) counters***
  integer  ::  i, j, k, n
  
  
  !!** this is where we read and write the data
  real, pointer, dimension(:,:,:,:)  :: solnData
  
  !!***vectors that store cell dimensions **
  real, allocatable, dimension(:) :: x, y, z
  
  !!***coordinates of grid cells***
  real :: xx,  yy, zz
  
  !!***variables that you set as you go 
  real   ::   vx, vy, vz, p, rho, metals
  real   ::   e, ek, gam, T
  
  !!*** sizes in each of the three direction
  integer :: sizeX,sizeY,sizeZ
  integer :: istat
  
  !! This says that we are grabing the guard cells when we grab the block
  logical :: gcell = .true.
  
  
  
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  
  !This is a part dedicated to the PrimordialChemistry. This is followed
  !from the Cellular problem and the unitTest problem
  real :: abar, zbar
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real :: temp_zone, rho_zone, vel_zone
  real ::  smallx
  real :: ptot, eint, etot
  real, dimension(EOS_NUM) :: eosData
  real :: hpA, hmA, dpA, hepA, heppA, h2pA, elecA
  
  
  !  print *,'Simulation InitBlock'
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData)
  
  !!***get coordinates of this block***
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(x(sizeX),stat=istat)
  allocate(y(sizeY),stat=istat)
  allocate(z(sizeZ),stat=istat)
  x = 0.0
  y = 0.0
  z = 0.0
  
  
  if (NDIM==3) call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,z,sizeZ)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,y,sizey)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,x,sizex)
  
  !Setup initial composition of all species
  !Question, Cellular had some small fraction for each species, can this not be zero?
  smallx = 1.0e-30;  !Trying something small
  massFraction(:) = smallx
  
! this initialization is basically free and arbitrarily adjustable,
! just make sure that the electron fraction is initialized properly to guarantee that the gas is neutral  !!!!!!
  if (HP_SPEC > 0) massFraction(HP_SPEC) = max(sim_xH*sim_xHP, smallx)
  if (HM_SPEC > 0) massFraction(HM_SPEC) = max(sim_xHM*massFraction(HP_SPEC)*sim_c_temp**0.88, smallx)
  if (HEP_SPEC > 0) massFraction(HEP_SPEC) = max(sim_xHE*4.*sim_xHeP, smallx)
  if (HEPP_SPEC > 0) massFraction(HEPP_SPEC) = max(sim_xHE*4.*sim_xHePP, smallx)

  if (H2P_SPEC > 0) massFraction(H2P_SPEC) = max(sim_xH2P*2.*massFraction(HP_SPEC)*sim_c_temp**1.8, smallx)
  if (H2_SPEC > 0) massFraction(H2_SPEC) = max(sim_xH2*2.*sim_xH*301.**5.1, smallx)

  if (H_SPEC > 0) massFraction(H_SPEC) = max(sim_xH-massFraction(HP_SPEC)-massFraction(HM_SPEC)- &
                                            & massFraction(H2P_SPEC)-massFraction(H2_SPEC),smallx)
  if (HE_SPEC > 0) massFraction(HE_SPEC) = max(sim_xHE-massFraction(HEP_SPEC)-massFraction(HEPP_SPEC), smallx)


  call Multispecies_getProperty(HP_SPEC,A,hpA)
  call Multispecies_getProperty(HM_SPEC,A,hmA)
  call Multispecies_getProperty(HEP_SPEC,A,hepA)
  call Multispecies_getProperty(HEPP_SPEC,A,heppA)
  call Multispecies_getProperty(H2P_SPEC,A,h2pA)
  call Multispecies_getProperty(ELEC_SPEC,A,elecA)

! guarantee that gas is neutral
  if (ELEC_SPEC > 0) massFraction(ELEC_SPEC) = max( elecA*(massFraction(HP_SPEC)/hpA - massFraction(HM_SPEC)/hmA &
 				        	& + massFraction(HEP_SPEC)/hepA + 2.*massFraction(HEPP_SPEC)/heppA + massFraction(H2P_SPEC)/h2pA), smallx)


  !***************************************************
  !!***Now loop over all cells in the current block***
  !  
  zz=0.
  do k = 1, sizeZ
     if(NDIM==3) zz = z(k)
     
     do j = 1, sizeY
        yy = y(j)
        
        do i = 1, sizeX
           xx = x(i)

           rho = sim_c_den*sim_contrast
           T = 5.
           
!           if(xx .le. 0.225*3.0856e18) then ! 1D case
           if(sqrt(xx**2+yy**2.+zz**2.) .le. 0.225*3.0856e18) then

             rho= sim_c_den 
             T= sim_c_temp 
           endif

           p=rho*T*sim_gasConst !!This will stay the same for pressure
           
           !Getting ready to find Gamma and some other Eos stuff
           eosData(EOS_TEMP) = T
           eosData(EOS_DENS) = rho
           eosData(EOS_PRES) = p
           
           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
           
           
           ptot = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)
           
           vx   = 0.
           vy   = 0.
           vz   = 0.
           ek   = 0.
           
           
           solnData(DENS_VAR,i,j,k) = eosData(EOS_DENS)
           solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
           solnData(VELX_VAR,i,j,k) = vx
           solnData(VELY_VAR,i,j,k) = vy
           solnData(VELZ_VAR,i,j,k) = vz
           solnData(GAME_VAR,i,j,k) = eosData(EOS_GAMC)
           solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
           solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
           solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT) + ek
           solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
           
           do n=SPECIES_BEGIN,SPECIES_END
              solnData(n, i,j,k) = massFraction(n)
           enddo

        enddo
     enddo
  enddo
  
  
  call Grid_releaseBlkPtr(blockID,solnData)
  
  
  deallocate(x)
  deallocate(y)
  deallocate(z)
  
  return
  
end subroutine Simulation_initBlock

