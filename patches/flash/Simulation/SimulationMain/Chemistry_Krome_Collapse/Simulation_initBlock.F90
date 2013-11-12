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
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData
#KROME_specs_block_vars
  
  
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
  massFraction(:) = smallx


#KROME_specs_block_if

#KROME_specs_block_prop

#KROME_specs_block_e


  !***************************************************
  !!***Now loop over all cells in the current block***
  !  
  xx=0.
  yy=0.
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

