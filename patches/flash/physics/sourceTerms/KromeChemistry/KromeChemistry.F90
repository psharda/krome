!!****f* source/physics/sourceTerms/KromeChemistry/KromeChemistry
!!
!! NAME
!!
!! KromeChemistry
!!
!! SYNOPSIS
!!
!!  KromeChemistry(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt)
!!
!!
!!
!! DESCRIPTION
!!  Apply a cooling operator on the list of blocks provided as input
!!
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!!
!!
!!
!!******


subroutine KromeChemistry (blockCount, blockList, dt)

  implicit none
  
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt

return
end subroutine KromeChemistry
