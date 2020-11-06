module krome_commons
  implicit none

#KROME_header
#KROME_species_index
#KROME_parameters

#KROME_atom_index

  !cooling index
#KROME_cool_index

  !heating index
#KROME_heat_index

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
#KROME_implicit_arr_r
#KROME_implicit_arr_p

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

#IFKROME_useDust
  !commons for dust
  real*8::krome_dust_partner_ratio(ndust),krome_dust_partner_ratio_inv(ndust)
  real*8::krome_dust_partner_mass(ndustTypes)
  real*8::krome_dust_asize(ndust),krome_dust_T(ndust)
  real*8::krome_dust_asize2(ndust),krome_dust_aspan(ndust)
  real*8::krome_dust_asize3(ndust),krome_grain_rho(ndustTypes)
  real*8::krome_dust_Tbind(ndust)
  real*8,allocatable::dust_Qabs(:,:),dust_Qabs_E(:)
  real*8,allocatable::dust_intBB(:,:),dust_intBB_sigma(:,:)
  real*8,allocatable::dust_intBB_dT(:,:)
  logical::Qabs_allocated(ndustTypes)
  real*8::dust_Qabs_interp(nPhotoBins,ndust)
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

  ! Draine dust absorption data loaded from file, via load_kabs
  ! in krome_photo module
  real*8::find_Av_draine_kabs(nPhotoBins)

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties
#KROME_opt_variables

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
#KROME_user_commons

#IFKROME_useSemenov
  !define the common Ebinding 
  real*8::Ebinding(nspec)
#ENDIFKROME

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
  real*8::dust2gas_ratio

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

#IFKROME_useCoolingZCIE
!data for metal cooling from table in the presence of UV
  integer,parameter::coolZCIEn1=81
  integer,parameter::coolZCIEn2=81
  integer,parameter::coolZCIEn3=81
  real*8::coolZCIEx1(coolZCIEn1),coolZCIEx2(coolZCIEn2),coolZCIEx3(coolZCIEn3)
  real*8::coolZCIEixd1(coolZCIEn1-1),coolZCIEixd2(coolZCIEn2-1)
  real*8::coolZCIEixd3(coolZCIEn3-1)
  real*8::coolZCIEy(coolZCIEn1,coolZCIEn2,coolZCIEn3)
  real*8::heatZCIEy(coolZCIEn1,coolZCIEn2,coolZCIEn3)
  real*8::coolZCIEx1min,coolZCIEx1max
  real*8::coolZCIEx2min,coolZCIEx2max
  real*8::coolZCIEx3min,coolZCIEx3max
  real*8::coolZCIEdvn1,coolZCIEdvn2,coolZCIEdvn3
#ENDIFKROME

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

#IFKROME_dust_table_2D
  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)
#ENDIFKROME

#IFKROME_dust_table_3D
  integer,parameter::dust_tab_kmax=10
  real*8::dust_mult_AvVariable,dust_tab_AvVariable(dust_tab_kmax)
  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax, dust_tab_kmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax, dust_tab_kmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax, dust_tab_kmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax, dust_tab_kmax)
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
  !$omp threadprivate(last_coe)

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

#IFKROME_useCoolingOH
  !data for OH cooling
  integer,parameter::coolOHn1=29
  integer,parameter::coolOHn2=21
  integer,parameter::coolOHn3=21
  real*8::coolOHx1(coolOHn1),coolOHx2(coolOHn2),coolOHx3(coolOHn3)
  real*8::coolOHixd1(coolOHn1-1),coolOHixd2(coolOHn2-1),coolOHixd3(coolOHn3-1)
  real*8::coolOHy(coolOHn3,coolOHn2,coolOHn1)
  real*8::coolOHx1min,coolOHx1max
  real*8::coolOHx2min,coolOHx2max
  real*8::coolOHx3min,coolOHx3max
  real*8::coolOHdvn1,coolOHdvn2,coolOHdvn3
#ENDIFKROME

#IFKROME_useCoolingH2O
  !data for H2O cooling
  integer,parameter::coolH2On1=29
  integer,parameter::coolH2On2=21
  integer,parameter::coolH2On3=21
  real*8::coolH2Ox1(coolH2On1),coolH2Ox2(coolH2On2),coolH2Ox3(coolH2On3)
  real*8::coolH2Oixd1(coolH2On1-1),coolH2Oixd2(coolH2On2-1),coolH2Oixd3(coolH2On3-1)
  real*8::coolH2Oy(coolH2On3,coolH2On2,coolH2On1)
  real*8::coolH2Ox1min,coolH2Ox1max
  real*8::coolH2Ox2min,coolH2Ox2max
  real*8::coolH2Ox3min,coolH2Ox3max
  real*8::coolH2Odvn1,coolH2Odvn2,coolH2Odvn3
#ENDIFKROME

#IFKROME_useCoolingHCN
  !data for HCN cooling
  integer,parameter::coolHCNn1=30
  integer,parameter::coolHCNn2=21
  integer,parameter::coolHCNn3=21
  real*8::coolHCNx1(coolHCNn1),coolHCNx2(coolHCNn2),coolHCNx3(coolHCNn3)
  real*8::coolHCNixd1(coolHCNn1-1),coolHCNixd2(coolHCNn2-1),coolHCNixd3(coolHCNn3-1)
  real*8::coolHCNy(coolHCNn3,coolHCNn2,coolHCNn1)
  real*8::coolHCNx1min,coolHCNx1max
  real*8::coolHCNx2min,coolHCNx2max
  real*8::coolHCNx3min,coolHCNx3max
  real*8::coolHCNdvn1,coolHCNdvn2,coolHCNdvn3
#ENDIFKROME

  !xsecs from file variables
#KROME_xsecs_from_file

  ! Gibbs free energy data from file variables
#KROME_GFE_from_file

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

#IFKROME_useCoolingGH
  ! Photo reaction rates relevant for Gnedin-Hollon cooling/heating function
  real*8::PLW,PHI,PHeI,PCVI
  !$omp threadprivate(PLW,PHI,PHeI,PCVI)
#ENDIFKROME

  !partition functions common vars
#KROME_var_parts

  !verbatim reactions
  character*50::reactionNames(nrea)

#IFKROME_hasStoreOnceRates
  !store evaluate once rate
  real*8::rateEvaluateOnce(nrea)
  !$omp threadprivate(rateEvaluateOnce)
#ENDIFKROME

end module krome_commons
