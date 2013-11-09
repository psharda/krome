!!****if* source/Simulation/SimulationMain/Chemistry_Krome_Collapse/Simulation_init
!!
!! NAME
!!  
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Krome collapse test
!!  
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()


!!***used modules from FLASH.***
  
  use Simulation_data
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get,RuntimeParameters_getPrev
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
#include "Flash.h"
#include "constants.h"




  real :: nblockx,nblocky,nblockz
  real :: sim_lrefine_max
 
!! These variables check to see if the Geometry is right 
  integer :: meshGeom
  logical :: validGeom


 

     print *,'Getting Physical Constants'
     call PhysicalConstants_get("ideal gas constant",sim_gasConst)
     
     call RuntimeParameters_get("xmin", sim_xMin)
     call RuntimeParameters_get("ymin", sim_yMin)
 
     call RuntimeParameters_get("xmax", sim_xMax)
     call RuntimeParameters_get("ymax", sim_yMax)

     print *,'Getting RuntimeParameters'

     call RuntimeParameters_get("sim_nblockx",nblockx)
     call RuntimeParameters_get("sim_nblocky",nblocky)

     if(NDIM == 3) then
        call RuntimeParameters_get("zmin", sim_zMin)
        call RuntimeParameters_get("zmax", sim_zMax)
        call RuntimeParameters_get("sim_nblockz", nblockz )
     endif


     call RuntimeParameters_get("sim_c_temp", sim_c_temp)
     call RuntimeParameters_get("sim_c_den", sim_c_den)
     call RuntimeParameters_get("sim_contrast", sim_contrast)

#KROME_specs_init
   
     call RuntimeParameters_get("smallx", smallx)

    return

end subroutine Simulation_init






