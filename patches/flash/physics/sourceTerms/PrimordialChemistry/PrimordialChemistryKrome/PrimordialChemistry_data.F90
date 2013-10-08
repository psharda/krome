!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryKrome/PrimordialChemistry_data
!!
!! NAME
!!
!!  PrimordialChemistry_data
!!
!!
!! SYNOPSIS
!!
!! use PrimordialChemistry_data
!!
!! DESCRIPTION
!!
!! Store the data for the Globular PrimordialChemistry unit
!!
!!
!!***

module PrimordialChemistry_data
    
  implicit none
  !! Put all the required data here
      
#include "constants.h"
#include "Flash.h"

  logical, save :: pchem_usePrimordialChemistry
  real, save    :: mp
  integer, save, dimension(NSPECIES) :: specieMap
  real, save, dimension(NSPECIES)    :: amu

end module PrimordialChemistry_data
