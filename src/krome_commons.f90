module krome_commons
  implicit none

#KROME_header
#KROME_species_index
#KROME_parameters

  real*8::arr_k(nrea)
  real*8::nold(nspec)

  !commons for rate tables
  integer,parameter::ktab_n=int(1e4)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx
  
  !commons for implicit RHS
  integer::arr_p1(nrea), arr_p2(nrea), arr_p3(nrea), arr_p4(nrea)
  integer::arr_r1(nrea), arr_r2(nrea), arr_r3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for dust
  real*8::krome_dust_partner_ratio(ndust), krome_dust_partner_ratio_inv(ndust)
  real*8::krome_dust_partner_idx(ndustTypes)
  real*8::krome_dust_partner_mass(ndustTypes)
  real*8::krome_dust_asize(ndust),krome_dust_T(ndust)
end module krome_commons
