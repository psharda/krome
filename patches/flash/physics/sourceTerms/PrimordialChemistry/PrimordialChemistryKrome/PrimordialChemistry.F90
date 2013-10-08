!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryKrome/PrimordialChemistry
!!
!! NAME
!!
!! PrimordialChemistry
!!
!! SYNOPSIS
!!
!!  call PrimordialChemistry( integer, intent(IN)  :: blockCount,
!!      integer(:), intent(IN) :: blockList,
!!      real, intent(IN)  :: dt)
!!
!! DESCRIPTION
!!
!! Apply chemistry to all blocks in specified list
!!
!! ARGUMENTS
!!
!! blockCount -- dimension of blockList
!! blockList  -- array of blocks which should receive chemistry
!! dt --     passed to the internal pchem_burner module
!!
!! PARAMETERS
!!
!! usePrimordialChemistry  -- Boolean, True. Turns on chemistry module
!!
!! NOTES
!!
!! Daniel Seifried 2013
!!
!!***

subroutine PrimordialChemistry(blockCount, blockList, dt)
   
  use PrimordialChemistry_data
  use Simulation_data
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getBlkPtr, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped, Eos


  use krome_main  

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  
  !args
  integer, INTENT(in)    :: blockCount
  integer, INTENT(in), DIMENSION(blockCount) :: blockList
  real, intent(IN)    :: dt
  
  !locals
  integer     :: i, j, k, n
  integer     :: blockID, thisBlock
  
  real, dimension(NSPECIES) :: x, xn
  real    :: tmp, rho
  
  integer, dimension(2,MDIM)  :: blkLimits, blkLimitsGC
  
  logical :: getGuardCells = .true.
 
  real, pointer, dimension(:,:,:,:) :: solnData
  
  ! -------------------- Check if chemistry is requested in runtime parameter
  if (.not. pchem_usePrimordialChemistry) return
  
  !--- Off to the races, 
  
  call Timers_start("chemistry")
  

  ! loop over list of blocks passed in
  do thisBlock = 1, blockCount
     
     blockID = blockList(thisBlock)
     ! get dimensions/limits and coordinates
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Get a pointer to solution data
     call Grid_getBlkPtr(blockID,solnData)

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              
              tmp = solnData(TEMP_VAR,i,j,k)
              rho = solnData(DENS_VAR,i,j,k)

              do n = 1, NSPECIES
                 x(n) = solnData(specieMap(n),i,j,k)
                 xn(n) = x(n) * rho / mp / amu(n)
              enddo
              
              !Do the chemistry
              call krome(xn,tmp,dt)

              !map the species back
              do n=1,NSPECIES
                 x(n) = xn(n) * mp * amu(n) / rho
                 solnData(specieMap(n),i,j,k) = x(n)
              enddo

              solnData(TEMP_VAR,i,j,k) = tmp
              
           enddo
        enddo
     enddo
     
     call Grid_releaseBlkPtr(blockID,solnData)

     call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockID)

  enddo
  
  call Timers_stop("chemistry")
  
  return
  
end subroutine PrimordialChemistry
