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

extern double krome_get_table_Tdust(double *x,double *Tgas);
extern double *krome_convert_xmoc(double** xmoc, int *imap);
extern void krome_return_xmoc(double *x, int *imap, double **xmoc);
extern double krome_num2col(double num, double *x, double Tgas);
extern void krome_print_phys_variables();
#IFKROME_useXrays
extern void krome_set_J21xray(double xarg);
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
extern double *krome_get_dust_distribution();
extern void krome_set_dust_distribution(double *arg);
extern double* krome_get_dust_size();
extern void krome_set_dust_size(double *arg);
extern void krome_set_Tdust(double arg);
extern void krome_set_Tdust_array(double *arr);
extern double krome_get_averaged_Tdust();
extern void krome_scale_dust_distribution(double xscale);
extern double *krome_get_Tdust();
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
extern void krome_set_photoBinJ(double *phbin);
extern void krome_set_photobinE_lr(double *phbinleft, double *phbinright);
extern void krome_set_photobinE_moc(double *binPos, double *binWidth);
extern void krome_set_photobinE_lin(double lower, double upper);
extern void krome_set_photobinE_log(double lower, double upper);
extern double *krome_get_photoBinJ();
extern double *krome_get_photoBinE_left();
extern double *krome_get_photoBinE_right();
extern double *krome_get_photoBinE_mid();
extern double *krome_get_photoBinE_delta();
extern double *krome_get_photoBin_rates();
extern double *krome_get_xsec(int idx);
extern double *krome_get_photoBin_heats();
extern void krome_photoBin_scale(double xscale);
extern void krome_photoBin_scale_array(double *xscale);
extern void krome_photoBin_restore();
extern void krome_photoBin_store();
// extern void krome_load_photoBin_file(char *fname);
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern void krome_set_photoBin_HMlog(double lower_in, double upper_in);
extern void krome_set_photoBin_BBlin(double lower, double upper, double Tbb);
extern void krome_set_photoBin_BBlog(double lower, double upper, double Tbb);
extern void krome_set_photoBin_BBlog_auto(double Tbb);
extern void krome_set_photoBin_draineLin(double lower, double upper);
extern void krome_set_photoBin_draineLog(double lower, double upper);
extern void krome_set_photoBin_J21lin(double lower, double upper);
extern void krome_set_photoBin_J21log(double lower, double upper);
extern double *krome_get_opacity(double *x, double Tgas);
extern double *krome_get_opacity_size(double *x, double Tgas, double csize);
extern void krome_dump_Jflux(int nfile);
#ENDIFKROME
extern double *krome_get_coef(double Tgas);
extern double krome_get_mu_x(double *xin);
extern double krome_get_gamma_x(double *xin, double inTgas);
extern void krome_consistent_x(double *x);
extern double *krome_n2x(double *n, double rhogas);
extern double *krome_x2n(double *x, double rhogas);
extern void krome_thermo(double *x, double *Tgas, double dt);
#IFKROME_use_heating
extern double krome_get_heating(double *x, double inTgas);
extern double *krome_get_heating_array(double *x, double inTgas);
#ENDIFKROME
#IFKROME_use_cooling
extern double krome_get_cooling(double *x, double inTgas);
extern double *krome_get_cooling_array(double *x, double inTgas);
extern void krome_plot_cooling(double *n);
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern void krome_dump_cooling(double *n, double Tgas, int nfile_in);
#ENDIFKROME
extern void krome_conserveLin_x(double *x, double *ref);
extern double *krome_conserveLinGetRef_x(double *x);
extern double *krome_conserve(double *x, double *xi);
extern double krome_get_gamma(double *x, double Tgas);
extern int *krome_get_zatoms();
extern double krome_get_mu(double *x);
// extern char **krome_get_rnames();
extern double *krome_get_mass();
extern double *krome_get_imass();
extern double krome_get_Hnuclei(double *x);
extern double *krome_get_charges();
// extern char **krome_get_names();
extern int krome_get_index(char *name);
extern double krome_get_rho(double *n);
extern void krome_scale_Z(double *n, double Z);
extern void krome_set_Z(double xarg);
extern void krome_set_clump(double xarg);
extern double krome_get_electrons(double *x);
extern void krome_print_best_flux(double *xin, double Tgas, int nbest);
extern void krome_print_best_flux_frac(double *xin, double Tgas, double frac);
extern void krome_print_best_flux_spec(double *xin, double Tgas, int nbest, int idx_find);
extern double *krome_get_flux(double *n, double Tgas);
extern void krome_explore_flux(double *x, double Tgas, int ifile, double xvar);
extern double *krome_get_qeff();
#IFKROME_useStars
/* Here is another example of a Fortran subroutine which accepts optional arguments. */
extern double *krome_stars_coe(double *x, double rho, double Tgas, double *y_in, double *zz_in);
extern double *krome_stars_energy(double *x, double rho, double Tgas, double *k);
#ENDIFKROME
extern void krome_dump_flux(double *n, double Tgas, int nfile);
extern void krome_dump_rates(double inTmin, double inTmax, int imax, int funit);
extern void krome_get_info(double *x, double Tgas);

#endif
