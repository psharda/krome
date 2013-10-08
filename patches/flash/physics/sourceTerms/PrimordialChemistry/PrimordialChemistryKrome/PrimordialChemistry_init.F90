!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryKrome/PrimordialChemistry_init
!!
!! NAME
!!
!!  PrimordialChemistry_init
!!
!! SYNOPSIS
!!
!!  PrimordialChemistry_init()
!!
!! 
!! DESCRIPTION
!!
!!  Initalizes various runtime paramters for PrimordialChemistry
!!
!!  ARGUMENTS
!!
!!  PARAMETERS
!!  
!!   usePrimordialChemistry -- Boolean, True. Turns on PrimordialChemistry module
!!
!!***

subroutine PrimordialChemistry_init()
  
   use PrimordialChemistry_data
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use PhysicalConstants_interface, ONLY: PhysicalConstants_get
   use Multispecies_interface, ONLY : Multispecies_getProperty

   use krome_user

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

   integer :: n
   character :: krome_species_names*16(NSPECIES)
   
   call RuntimeParameters_get("usePrimordialChemistry", pchem_usePrimordialChemistry)
   call PhysicalConstants_get("proton mass",mp)

   krome_species_names = krome_get_names()   

   do n = 1, NSPECIES
     call pchem_mapNetworkToSpecies(krome_species_names(n),specieMap(n))
     call Multispecies_getProperty(specieMap(n), A, amu(n))
   enddo

   return

end subroutine PrimordialChemistry_init
