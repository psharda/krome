module krome_commons
  implicit none

#KROME_header
#KROME_species_index
#KROME_parameters

  !cooling index
#KROME_cool_index

  !heating index
#KROME_heat_index

  real*8::arr_k(nrea)
  real*8::jac_nold(nspec),jac_dnold(nspec),jac_dn(nspec)
  !$omp threadprivate(jac_nold,jac_dnold,jac_dn)
    
  !commons for rate tables
  !modify ktab_n according to the requested precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)
  
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
  real*8,allocatable::dust_Qabs(:,:),dust_Qabs_E(:),dust_intBB(:,:)
  real*8::xdust(ndust)
  integer::krome_dust_partner_idx(ndustTypes),dust_Qabs_nE
  integer,parameter::dust_nT=int(1e3)
  real*8::dust_intBB_Tbb(dust_nT),dust_cooling
  
  !commons for frequency bins
#KROME_photobins_array

  !commons for dust optical properties
#KROME_opt_variables

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
#KROME_user_commons

  !commons for anytab
#KROME_vars_anytab

  !commons for H2_opacity
  real*8::arrH2esc_Tgas(13),arrH2esc_ntot(10),arrH2esc(10,13)
  real*8::xmulH2esc,ymulH2esc

  !physical commons
#KROME_phys_commons

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !commons for Mayer opacity table
  integer,parameter::mayern=15,mayerm=23
  real*8::mayer_x(mayern),mayer_y(mayerm)
  real*8::mayer_z(mayern,mayerm),mayer_xmul,mayer_ymul
  
end module krome_commons
