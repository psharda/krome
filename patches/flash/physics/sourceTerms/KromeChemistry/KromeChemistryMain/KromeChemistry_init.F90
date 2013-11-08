!!****if* source/physics/sourceTerms/KromeChemistry/KromeChemistry_init
!!
!! NAME
!!
!!  KromeChemistry_init
!!
!! SYNOPSIS
!!
!!  KromeChemistry_init()
!!
!! 
!! DESCRIPTION
!!
!!  Initalizes various runtime paramters for KromeChemistry
!!
!!  ARGUMENTS
!!
!!  PARAMETERS
!!  
!!   useKromeChemistry -- Boolean, True. Turns on KromeChemistry module
!!
!! D. Seifried and KROME team 2013
!!
!!***

subroutine KromeChemistry_init()
  
   use KromeChemistry_data
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use PhysicalConstants_interface, ONLY: PhysicalConstants_get
   use Multispecies_interface, ONLY : Multispecies_getProperty

   use krome_user
   use krome_main

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

   integer :: n
   character :: krome_species_names*16(NSPECIES)
   
   call RuntimeParameters_get("useKromeChemistry", pchem_useKromeChemistry)
   call PhysicalConstants_get("proton mass",mp)

   call krome_init()

   krome_species_names = krome_get_names()   

   do n = 1, NSPECIES
     call pchem_mapNetworkToSpecies(krome_species_names(n),specieMap(n))
     call Multispecies_getProperty(specieMap(n), A, amu(n))
   enddo

   return

end subroutine KromeChemistry_init
