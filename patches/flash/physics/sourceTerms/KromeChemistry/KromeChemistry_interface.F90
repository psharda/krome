!!****h* source/physics/sourceTerms/KromeChemistry/KromeChemistry_interface
!! NAME
!!
!! KromeChemistry_interface
!!
!! SYNOPSIS
!!
!! use KromeChemistry_interface
!!
!! This is the header file for the kromechemistry module that defines its
!! public interfaces.
!!***

Module KromeChemistry_interface
#include "constants.h"
#include "Flash.h"

  interface 
     subroutine KromeChemistry_computeDt (blockID,  blkLimits, blkLimitsGC, solnData, dt_check, dt_minloc)
       implicit none
       integer, intent(IN) :: blockID
       integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
       real,INTENT(INOUT)  :: dt_check
       integer,INTENT(INOUT) :: dt_minloc(5)
       real, pointer :: solnData(:,:,:,:)
     end subroutine KromeChemistry_computeDt
  end interface
  
  interface 
     subroutine KromeChemistry(blockCount,blockList,dt)
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount), intent(IN) :: blockList
       real,intent(IN) :: dt
     end subroutine KromeChemistry
  end interface
  
  interface 
     subroutine KromeChemistry_finalize()
     end subroutine KromeChemistry_finalize
  end interface
  
  interface 
     subroutine KromeChemistry_init()
     end subroutine KromeChemistry_init
  end interface
  
end Module KromeChemistry_interface
