module krome_commons
  implicit none

#KROME_header
#KROME_species_index
#KROME_parameters

  real*8::arr_k(nrea)
  real*8::jac_nold(nspec),jac_dnold(nspec),jac_dn(nspec)
    
  !commons for rate tables
  !modify ktab_n according to the requested precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx
  
  !commons for implicit RHS
#KROME_implicit_arr_r
#KROME_implicit_arr_p

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for dust
  real*8::krome_dust_partner_ratio(ndust), krome_dust_partner_ratio_inv(ndust)
  real*8::krome_dust_partner_mass(ndustTypes)
  real*8::krome_dust_asize(ndust),krome_dust_T(ndust)
  real*8::krome_dust_asize2(ndust),krome_dust_aspan(ndust)
  real*8::krome_dust_asize3(ndust),krome_grain_rho
  integer::krome_dust_partner_idx(ndustTypes)

  !commons for photoionization
#KROME_photo_variables
#KROME_photoheating_variables

  !commons for dust optical properties
#KROME_opt_variables

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0

end module krome_commons
