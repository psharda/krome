module krome_constants
  implicit none

  real*8,parameter:: boltzmann_eV   = 8.617332478d-5 ! eV / K 
  real*8,parameter:: boltzmann_J    = 1.380648d-23   ! J / K    
  real*8,parameter:: boltzmann_erg  = 1.380648d-16   ! erg / K     

  real*8,parameter:: planck_eV   =  4.135667516d-15! eV s 
  real*8,parameter:: planck_J    = 6.62606957d-34   ! J s    
  real*8,parameter:: planck_erg  = 6.62606957d-27   ! erg s    

  real*8,parameter:: gravity  = 6.674d-8 ! cm3 / g / s2       
  real*8,parameter:: e_mass   = 9.10938188d-28  ! g
  real*8,parameter:: p_mass   = 1.67262158d-24  ! g
  real*8,parameter:: clight = 2.99792458e10 !cm/s

  real*8,parameter:: pi = 3.14159 !PI
  
  real*8,parameter:: eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter:: seconds_per_year = 365d0*24d0*3600d0 !yr -> s

  real*8,parameter:: kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass
  
  real*8,parameter:: pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter:: exp_planck =  planck_erg / boltzmann_erg !s*K
  real*8,parameter:: stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4

  real*8,parameter:: N_avogadro = 6.0221d23 !#
  real*8,parameter:: Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter:: Rgas_kJ = 8.3144621d-3 !kJ/K/mol


  
end module krome_constants
