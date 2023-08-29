/* C header for krome_user.f90/the krome_user module.
  Variables for which KROME will return an updated value need to be passed
  by reference using an argument pointer, e.g., "*x".
  Arrays as input should also be passed by reference/argument pointer.
  Passing 2-D arrays as arguments use an array of pointers to an array (i.e., "**x").
  To return an array from a function, it must return an argument pointer,
  e.g., "extern double *functionName()".
*/
#ifndef KROME_USER_H_
#define KROME_USER_H_

#KROME_species

#KROME_cool_index

#KROME_heat_index

#KROME_common_alias

#KROME_constant_list

#KROME_user_commons_functions

#KROME_set_get_phys_functions

#KROME_cooling_functions

extern double krome_get_table_tdust(double *x,double *Tgas);
extern double krome_num2col(double num, double *x, double Tgas);
extern void krome_print_phys_variables();
extern void krome_set_mpi_rank(int rank);
#IFKROME_useXrays
extern void krome_set_j21xray(double xarg);
#ENDIFKROME
#IFKROME_use_coolingZ
extern void krome_popcool_dump(double xvar, int nfile);
#ENDIFKROME
#IFKROME_useTabsTdust
extern double krome_get_Tdust(double *x, double Tgas);
#ENDIFKROME
#IFKROME_useDust
/* NOTE: The following Fortran void accepts optional arguments, which don't
   seem to be intrinsically possible in C. Let's make them required arguments
   for now.
*/
extern void krome_init_dust_distribution(double *x, double dust_gas_ratio,
    double alow_arg, double aup_arg, double phi_arg);
extern void krome_get_dust_distribution(double* krome_get_dust_distribution_var);
extern void krome_set_dust_distribution(double *arg);
extern void krome_get_dust_size(double* krome_get_dust_size_var);
extern void krome_set_dust_size(double *arg);
extern void krome_set_tdust(double arg);
extern void krome_set_tdust_array(double *arr);
extern double krome_get_averaged_Tdust();
extern void krome_scale_dust_distribution(double xscale);
extern void krome_get_tdust(double* krome_get_tdust_var);
extern void krome_set_surface(double *x, double xarg, int idx_base);
extern void krome_set_surface_norm(double *x, double xarg, int idx_base);
extern void krome_set_surface_array(double *x, double *xarr, int idx_base);
extern double krome_get_surface(double *x, int idx_base);
#ENDIFKROME
extern void krome_store(double *x, double Tgas, double dt);
extern void krome_restore(double *x, double *Tgas, double *dt);
extern void krome_thermo_on();
extern void krome_thermo_off();
#IFKROME_usePhotoBins
extern void krome_calc_photobins();
extern void krome_set_photobinj(double *phbin);
extern void krome_set_photobine_lr(double *phbinleft, double *phbinright, double tgas);
extern void krome_set_photobine_moc(double *binPos, double *binWidth, double tgas);
extern void krome_set_photobine_lin(double lower, double upper, double tgas);
extern void krome_set_photobine_log(double lower, double upper, double tgas);
extern void krome_get_photobinj(double* krome_get_photobinj_var);
extern void krome_get_photobine_left(double* krome_get_photobine_left_var);
extern void krome_get_photobine_right(double* krome_get_photobine_right_var);
extern void krome_get_photobine_mid(double* krome_get_photobine_mid_var);
extern void krome_get_photobine_delta(double* krome_get_photobine_delta_var);
extern void krome_get_photobin_rates(double* krome_get_photobin_rates_var);
extern void krome_get_xsec(int* idx,double* krome_get_xsec_var);
extern void krome_get_photobin_heats(double* krome_get_photobin_heats_var);
extern void krome_photobin_scale(double xscale);
extern void krome_photobin_scale_array(double *xscale);
extern void krome_photobin_restore();
extern void krome_photobin_store();
// extern void krome_load_photoBin_file(char *fname);
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern void krome_set_photobin_hmlog(double lower_in, double upper_in);
extern void krome_set_photobin_hmcustom(double lower_in,double upper_in,int* additive);
extern void krome_set_photobin_bblin(double lower, double upper, double Tbb);
extern void krome_set_photobin_bblog(double lower, double upper, double Tbb);
extern void krome_set_photobin_bblog_auto(double Tbb);
extern void krome_set_photobin_drainelin(double lower, double upper);
extern void krome_set_photobin_drainelog(double lower, double upper);
extern void krome_set_photobin_j21lin(double lower, double upper);
extern void krome_set_photobin_j21log(double lower, double upper);
extern double krome_get_photointensity(double energy);
extern void krome_get_opacity(double *x, double Tgas, double* krome_get_opacity_var);
extern void krome_get_opacity_size(double *x, double Tgas, double csize, double* krome_get_opacity_size_var);
extern void krome_get_opacity_size_d2g(double *x, double Tgas, double csize, double d2g, double* krome_get_opacity_size_d2g_var);
extern void krome_load_opacity_table();
extern void krome_dump_jflux(int nfile);
#ENDIFKROME
extern void krome_get_coef(double Tgas, double* x, double* krome_get_coef_var);
extern double krome_get_mu_x(double *xin);
extern double krome_get_gamma_x(double *xin, double inTgas);
extern void krome_consistent_x(double *x);
extern void krome_n2x(double *n, double rhogas, double* krome_n2x_var);
extern void krome_x2n(double *x, double rhogas, double* krome_x2n_var);
extern void krome_thermo(double *x, double *Tgas, double dt);
#IFKROME_use_heating
extern double krome_get_heating(double *x, double inTgas);
extern void krome_get_heating_array(double *x, double inTgas, double* krome_get_heating_array_var);
#ENDIFKROME
#IFKROME_use_cooling
extern double krome_get_cooling(double *x, double inTgas);
extern void krome_get_cooling_array(double *x, double inTgas, double* krome_get_cooling_array_var);
extern void krome_plot_cooling(double *n);
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern void krome_dump_cooling(double *n, double Tgas, int nfile_in);
#ENDIFKROME
extern void krome_conservelin_x(double *x, double *ref);
extern void krome_conservelingetref_x(double *x, double* krome_conserveLinGetRef_var);
extern void krome_conserve(double *x, double *xi, double* krome_conserve_var);
extern double krome_get_gamma(double *x, double Tgas);
extern void krome_get_zatoms(int* krome_get_zatoms_var);
extern double krome_get_mu(double *x);
// extern char **krome_get_rnames();
extern void krome_get_mass(double* krome_get_mass_var);
extern void krome_get_imass(double* krome_get_imass_var);
extern double krome_get_hnuclei(double *x);
extern void krome_get_charges(double* krome_get_charges_var);
// extern char **krome_get_names();
extern int krome_get_index(char *name);
extern double krome_get_rho(double *n);
extern void krome_scale_z(double *n, double Z);
extern void krome_set_z(double xarg);
extern void krome_set_dust_to_gas(double xarg);
extern void krome_set_clump(double xarg);
extern double krome_get_electrons(double *x);
extern void krome_print_best_flux(double *xin, double Tgas, int nbest);
extern void krome_print_best_flux_frac(double *xin, double Tgas, double frac);
extern void krome_print_best_flux_spec(double *xin, double Tgas, int nbest, int idx_find);
extern void krome_get_flux(double *n, double Tgas, double* krome_get_flux_var);
extern void krome_explore_flux(double *x, double Tgas, int ifile, double xvar);
extern void krome_get_qeff(double* krome_get_qeff_var);
#IFKROME_useStars
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern void krome_stars_coe(double *x, double rho, double Tgas, double *y_in, double *zz_in, double* krome_stars_coe_var);
extern void krome_stars_energy(double *x, double rho, double Tgas, double *k, double* krome_stars_energy_var);
#ENDIFKROME
extern void krome_dump_flux(double *n, double Tgas, int nfile);
extern void krome_dump_rates(double inTmin, double inTmax, int imax, int funit);
extern void krome_get_info(double *x, double Tgas);

#endif
