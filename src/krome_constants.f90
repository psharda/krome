module krome_constants
  implicit none

  double precision, parameter :: boltzmann_eV   = 8.617332478d-5 ! eV / K 
  double precision, parameter :: boltzmann_J    = 1.380648d-23   ! J / K    
  double precision, parameter :: boltzmann_erg  = 1.380648d-16   ! erg / K     

  double precision, parameter :: planck_eV   =  4.135667516d-15! eV s 
  double precision, parameter :: planck_J    = 6.62606957d-34   ! J s    
  double precision, parameter :: planck_erg  = 6.62606957d-27   ! erg / K     

  double precision, parameter :: gravity  = 6.674d-8 ! cm3 / g / s2       
  double precision, parameter :: e_mass   = 9.10938188d-28  ! g
  double precision, parameter :: p_mass   = 1.67262158d-24  ! g

  double precision, parameter :: pi = 3.14159 !PI
  
  double precision, parameter :: eV_to_erg = 1.60217646d-12 !eV -> erg
  double precision, parameter :: seconds_per_year = 365d0*24d0*3600d0 


end module krome_constants
