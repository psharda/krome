!!****if* source/physics/sourceTerms/KromeChemistry/KromeChemistry_data
!!
!! NAME
!!
!!  KromeChemistry_data
!!
!!
!! SYNOPSIS
!!
!! use KRomeChemistry_data
!!
!! DESCRIPTION
!!
!! Store the data for the Globular KromeChemistry unit
!!
!! D. Seifried and KROME Team 2013
!!
!!***

module KromeChemistry_data
    
  implicit none
  !! Put all the required data here
      
#include "constants.h"
#include "Flash.h"

  logical, save :: pchem_useKromeChemistry
  real, save    :: mp
  integer, save, dimension(NSPECIES) :: specieMap
  real, save, dimension(NSPECIES)    :: amu

end module KromeChemistry_data
