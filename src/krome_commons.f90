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

#IFKROME_useDust
  !commons for dust
  real*8::krome_dust_partner_ratio(ndust), krome_dust_partner_ratio_inv(ndust)
  real*8::krome_dust_partner_mass(ndustTypes)
  real*8::krome_dust_asize(ndust),krome_dust_T(ndust)
  real*8::krome_dust_asize2(ndust),krome_dust_aspan(ndust)
  real*8::krome_dust_asize3(ndust),krome_grain_rho(ndustTypes)
  real*8::krome_dust_Tevap(ndust)
  real*8,allocatable::dust_Qabs(:,:),dust_Qabs_E(:)
  real*8,allocatable::dust_intBB(:,:),dust_intBB_sigma(:,:)
  real*8,allocatable::dust_intBB_dT(:,:)
  real*8::xdust(ndust)
  integer::krome_dust_partner_idx(ndustTypes),dust_Qabs_nE
  integer,parameter::dust_nT=int(1e5)
  real*8::dust_intBB_Tbb(dust_nT),dust_cooling
  !logarithm of the maximum BB temperature in integral tables
  real*8,parameter::TbbMax=1d4,TbbMin=1d0
  real*8,parameter::TbbMult=(dust_nT-1)/(Tbbmax-Tbbmin)
#ENDIFKROME

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

#IFKROME_useOmukaiOpacity
  !commons for H2_opacity
  real*8::arrH2esc_Tgas(20),arrH2esc_ntot(10),arrH2esc(10,20)
  real*8::xmulH2esc,ymulH2esc
#ENDIFKROME

  !physical commons
#KROME_phys_commons

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar 
  real*8::total_Z

#IFKROME_useMayerOpacity
  !commons for Mayer opacity table
  integer,parameter::mayern=15,mayerm=23
  real*8::mayer_x(mayern),mayer_y(mayerm)
  real*8::mayer_z(mayern,mayerm),mayer_xmul,mayer_ymul
#ENDIFKROME

#IFKROME_useCoolingZCIENOUV
  integer,parameter::CoolZNOUVn=131,CoolZNOUVm=162
  real*8::CoolZNOUV_x(CoolZNOUVn),CoolZNOUV_y(CoolZNOUVm)
  real*8::CoolZNOUV_z(CoolZNOUVn,CoolZNOUVm),CoolZNOUV_xmul,CoolZNOUV_ymul
#ENDIFKROME

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
       / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

#IFKROME_useChemisorption
  !commons for chemisorption rate fit
  real*8,allocatable::dust_rateChem_PC(:),dust_rateChem_CP(:)
  real*8,allocatable::dust_rateChem_CC(:),dust_rateChem_x(:)
  real*8::dust_rateChem_xmin,dust_rateChem_dx
  real*8::dust_rateChem_xfact,dust_rateChem_invdx
  integer::dust_rateChem_xsteps
#ENDIFKROME

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)

#IFKROME_usePreDustExp
  !store the exponentials at constant Tdust
  real*8::dust_Ebareice_exp(2*nspec)
  real*8::dust_Ebareice23_exp(2*nspec)
#ENDIFKROME

#IFKROME_useCoolingCO
  !data for CO cooling
  integer,parameter::coolCOn1=40
  integer,parameter::coolCOn2=40
  integer,parameter::coolCOn3=40
  real*8::coolCOx1(coolCOn1),coolCOx2(coolCOn2),coolCOx3(coolCOn3)
  real*8::coolCOixd1(coolCOn1-1),coolCOixd2(coolCOn2-1),coolCOixd3(coolCOn3-1)
  real*8::coolCOy(coolCOn1,coolCOn2,coolCOn3)
  real*8::coolCOx1min,coolCOx1max
  real*8::coolCOx2min,coolCOx2max
  real*8::coolCOx3min,coolCOx3max
  real*8::coolCOdvn1,coolCOdvn2,coolCOdvn3
#ENDIFKROME

#IFKROME_useCoolingZCIE
  !data for metal CIE cooling
  integer,parameter::coolZCIEn1=81 
  integer,parameter::coolZCIEn2=81
  integer,parameter::coolZCIEn3=81
  real*8::coolZCIEx1(coolZCIEn1),coolZCIEx2(coolZCIEn2),coolZCIEx3(coolZCIEn3)
  real*8::coolZCIEixd1(coolZCIEn1-1),coolZCIEixd2(coolZCIEn2-1),coolZCIEixd3(coolZCIEn3-1)
  real*8::coolZCIEy(coolZCIEn3,coolZCIEn2,coolZCIEn1)
  real*8::heatZCIEy(coolZCIEn3,coolZCIEn2,coolZCIEn1)
  real*8::coolZCIEx1min,coolZCIEx1max
  real*8::coolZCIEx2min,coolZCIEx2max
  real*8::coolZCIEx3min,coolZCIEx3max
  real*8::coolZCIEdvn1,coolZCIEdvn2,coolZCIEdvn3
#ENDIFKROME

  !xsecs from file variables
#KROME_xsecs_from_file

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !partition functions common vars
#KROME_var_parts

end module krome_commons
