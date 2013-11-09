!!****if* source/Simulation/SimulationMain/Chemistry_Krome_Collapse/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for Krome collapse test problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! These are the physical constants we use
  real, save :: sim_gasConst

 real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
 real, save    :: sim_c_den, sim_c_temp
 real, save    :: sim_contrast

!!Fluid

#KROME_specs_data

 real, save :: smallx

end module Simulation_data
