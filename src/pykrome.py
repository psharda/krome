# -*- coding: utf-8 -*-
import numpy as np
import ctypes
import numpy.ctypeslib as npctypes

# define aliases for complicated variable types
int_byref = ctypes.POINTER(ctypes.c_int)
dble_byref = ctypes.POINTER(ctypes.c_double)
array_1d_int = npctypes.ndpointer(dtype=np.int,ndim=1,flags='CONTIGUOUS')
array_1d_double = npctypes.ndpointer(dtype=np.double,ndim=1,flags='CONTIGUOUS')
array_2d_double = npctypes.ndpointer(dtype=np.double,ndim=2,flags='CONTIGUOUS')

class PyKROME(object):
	# NOTE: Fortran is, by default, case-insensitive. As a side effect, all of the
	# symbols in `libkrome.so` are entirely lower-case (regardless of their form in
	# the Fortran files).

	def __init__(self, path_to_lib='.'):
		fortran = npctypes.load_library('libkrome.so',path_to_lib)

		# alias the Fortran functions and subroutines to `lib`.
		self.lib = fortran

		#KROME_species

		#KROME_cool_index

		#KROME_heat_index

		#KROME_common_alias

		#KROME_constant_list

		#KROME_user_commons_functions

		#KROME_set_get_phys_functions

		#KROME_cooling_functions

		# krome.f90 argument and return (result) types
		fortran.krome.restype = None
		fortran.krome_equilibrium.restype = None
		#IFKROME_useX
		fortran.krome.argtypes = [array_1d_double, dble_byref, dble_byref, dble_byref]
		fortran.krome_equilibrium.argtypes = [array_1d_double, dble_byref, dble_byref]
		#ELSEKROME
		fortran.krome.argtypes = [array_1d_double, dble_byref, dble_byref]
		fortran.krome_equilibrium.argtypes = [array_1d_double, dble_byref]
		#ENDIFKROME
		#IFKROME_ierr
		fortran.krome.argtypes.append(int_byref)
		#ENDIFKROME
		fortran.krome_init.restype = None
		fortran.krome_init.argtypes = None
		fortran.krome_get_coe.restype = array_1d_double
		fortran.krome_get_coe.argtypes = [array_1d_double, dble_byref]
		fortran.krome_get_coet.restype = array_1d_double
		fortran.krome_get_coet.argtypes = [dble_byref]
	
		#krome_user.f90 argument and return (result) types
		fortran.krome_get_table_tdust.restype = None
		fortran.krome_get_table_tdust.argtypes = [array_1d_double, dble_byref]
		fortran.krome_num2col.restype = dble_byref
		fortran.krome_num2col.argtypes = [ctypes.c_double, array_1d_double, ctypes.c_double]
		fortran.krome_print_phys_variables.restype = None
		fortran.krome_print_phys_variables.argtypes = None
		#IFKROME_useXrays
		fortran.krome_set_j21xray.restype = None
		fortran.krome_set_j21xray.argtypes = [ctypes.c_double]
		#ENDIFKROME
		#IFKROME_use_coolingZ
		fortran.krome_popcool_dump.restype = None
		fortran.krome_popcool_dump.argtypes = [ctypes.c_double, ctypes.c_int]
		#ENDIFKROME
		#IFKROME_useTabsTdust
		fortran.krome_get_tdust.restype = dble_byref
		fortran.krome_get_tdust.argtypes = [array_1d_double, ctypes.c_double]
		#ENDIFKROME
		#IFKROME_useDust
		# NOTE: The following Fortran subroutine accepts optional arguments, which don't
		# seem to be intrinsically possible in C (but what about Python?). Let's make
		# them required arguments for now.
		fortran.krome_init_dust_distribution.restype = None
		fortran.krome_init_dust_distribution.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_get_dust_distribution.restype = array_1d_double
		fortran.krome_get_dust_distribution.argtypes = None
		fortran.krome_set_dust_distribution.restype = None
		fortran.krome_set_dust_distribution.argtypes = [array_1d_double]
		fortran.krome_get_dust_size.restype = array_1d_double
		fortran.krome_get_dust_size.argtypes = None
		fortran.krome_set_dust_size.restype = None
		fortran.krome_set_dust_size.argtypes = [array_1d_double]
		fortran.krome_set_tdust.restype = None
		fortran.krome_set_tdust.argtypes = [ctypes.c_double]
		fortran.krome_set_tdust_array.restype = None
		fortran.krome_set_tdust_array.argtypes = [array_1d_double]
		fortran.krome_get_averaged_tdust.restype = ctypes.c_double
		fortran.krome_get_averaged_tdust.argtypes = None
		fortran.krome_scale_dust_distribution.restype = None
		fortran.krome_scale_dust_distribution.argtypes = [ctypes.c_double]
		fortran.krome_get_tdust.restype = array_1d_double
		fortran.krome_get_tdust.argtypes = None
		fortran.krome_set_surface.restype = None
		fortran.krome_set_surface.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int]
		fortran.krome_set_surface_norm.restype = None
		fortran.krome_set_surface_norm.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int]
		fortran.krome_set_surface_array.restype = None
		fortran.krome_set_surface_array.argtypes = [array_1d_double, array_1d_double, ctypes.c_int]
		fortran.krome_get_surface.restype = ctypes.c_double
		fortran.krome_get_surface.argtypes = [array_1d_double, ctypes.c_int]
		#ENDIFKROME
		fortran.krome_store.restype = None
		fortran.krome_store.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_restore.restype = None
		fortran.krome_restore.argtypes = [array_1d_double, dble_byref, dble_byref]
		fortran.krome_thermo_on.restype = None
		fortran.krome_thermo_on.argtypes = None
		fortran.krome_thermo_off.restype = None
		fortran.krome_thermo_off.argtypes = None
		#IFKROME_usePhotoBins
		fortran.krome_calc_photobins.restype = None
		fortran.krome_calc_photobins.argtypes = None
		fortran.krome_set_photobinj.restype = None
		fortran.krome_set_photobinj.argtypes = [array_1d_double]
		fortran.krome_set_photobine_lr.restype = None
		fortran.krome_set_photobine_lr.argtypes = [array_1d_double, array_1d_double]
		fortran.krome_set_photobine_moc.restype = None
		fortran.krome_set_photobine_moc.argtypes = [array_1d_double, array_1d_double]
		fortran.krome_set_photobine_lin.restype = None
		fortran.krome_set_photobine_lin.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobine_log.restype = None
		fortran.krome_set_photobine_log.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_get_photobinj.restype = array_1d_double
		fortran.krome_get_photobinj.argtypes = None
		fortran.krome_get_photobine_left.restype = array_1d_double
		fortran.krome_get_photobine_left.argtypes = None
		fortran.krome_get_photobine_right.restype = array_1d_double
		fortran.krome_get_photobine_right.argtypes = None
		fortran.krome_get_photobine_mid.restype = array_1d_double
		fortran.krome_get_photobine_mid.argtypes = None
		fortran.krome_get_photobine_delta.restype = array_1d_double
		fortran.krome_get_photobine_delta.argtypes = None
		fortran.krome_get_photobin_rates.restype = array_1d_double
		fortran.krome_get_photobin_rates.argtypes = None
		fortran.krome_get_xsec.restype = None
		fortran.krome_get_xsec.argtypes = [ctypes.c_int]
		fortran.krome_get_photobin_heats.restype = array_1d_double
		fortran.krome_get_photobin_heats.argtypes = None
		fortran.krome_photobin_scale.restype = None
		fortran.krome_photobin_scale.argtypes = [ctypes.c_double]
		fortran.krome_photobin_scale_array.restype = None
		fortran.krome_photobin_scale_array.argtypes = [array_1d_double]
		fortran.krome_photobin_restore.restype = None
		fortran.krome_photobin_restore.argtypes = None
		fortran.krome_photobin_restore.restype = None
		fortran.krome_photobin_restore.argtypes = None
		fortran.krome_photobin_store.restype = None
		fortran.krome_photobin_store.argtypes = None
		## extern void krome_load_photoBin_file(char *fname);
		# Here is another example of a Fortran subroutine which accepts optional arguments:
		fortran.krome_set_photobin_hmlog.restype = None
		fortran.krome_set_photobin_hmlog.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_bblin.restype = None
		fortran.krome_set_photobin_bblin.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_bblog.restype = None
		fortran.krome_set_photobin_bblog.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_bblog_auto.restype = None
		fortran.krome_set_photobin_bblog_auto.argtypes = [ctypes.c_double]
		fortran.krome_set_photobin_drainelin.restype = None
		fortran.krome_set_photobin_drainelin.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_drainelog.restype = None
		fortran.krome_set_photobin_drainelog.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_j21lin.restype = None
		fortran.krome_set_photobin_j21lin.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_set_photobin_j21log.restype = None
		fortran.krome_set_photobin_j21log.argtypes = [ctypes.c_double, ctypes.c_double]
		fortran.krome_get_opacity.restype = array_1d_double
		fortran.krome_get_opacity.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_get_opacity_size.restype = array_1d_double
		fortran.krome_get_opacity_size.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_dump_jflux.restype = None
		fortran.krome_dump_jflux.argtypes = [ctypes.c_int]
		#ENDIFKROME
		fortran.krome_get_coef.restype = array_1d_double
		fortran.krome_get_coef.argtypes = [ctypes.c_double]
		fortran.krome_get_mu_x.restype = ctypes.c_double
		fortran.krome_get_mu_x.argtypes = [array_1d_double]
		fortran.krome_get_gamma_x.restype = ctypes.c_double
		fortran.krome_get_gamma_x.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_consistent_x.restype = None
		fortran.krome_consistent_x.argtypes = [array_1d_double]
		fortran.krome_n2x.restype = array_1d_double
		fortran.krome_n2x.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_x2n.restype = array_1d_double
		fortran.krome_x2n.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_thermo.restype = None
		fortran.krome_thermo.argtypes = [array_1d_double, array_1d_double, ctypes.c_double]
		#IFKROME_use_heating
		fortran.krome_get_heating.restype = ctypes.c_double
		fortran.krome_get_heating.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_get_heating_array.restype = array_1d_double
		fortran.krome_get_heating_array.argtypes = [array_1d_double, ctypes.c_double]
		#ENDIFKROME
		#IFKROME_use_cooling
		fortran.krome_get_cooling.restype = ctypes.c_double
		fortran.krome_get_cooling.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_get_cooling_array.restype = array_1d_double
		fortran.krome_get_cooling_array.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_plot_cooling.restype = None
		fortran.krome_plot_cooling.argtypes = [array_1d_double]
		# Here is another example of a Fortran subroutine which accepts optional arguments
		# (and which are subsequently forced to be mandatory...for now).
		fortran.krome_dump_cooling.restype = None
		fortran.krome_dump_cooling.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int]
		#ENDIFKROME
		fortran.krome_conservelin_x.restype = None
		fortran.krome_conservelin_x.argtypes = [array_1d_double, array_1d_double]
		fortran.krome_conservelingetref_x.restype = ctypes.c_double
		fortran.krome_conservelingetref_x.argtypes = [array_1d_double]
		fortran.krome_conserve.restype = array_1d_double
		fortran.krome_conserve.argtypes = [array_1d_double, array_1d_double]
		fortran.krome_get_gamma.restype = ctypes.c_double
		fortran.krome_get_gamma.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_get_zatoms.restype = ctypes.c_int
		fortran.krome_get_zatoms.argtypes = None
		fortran.krome_get_mu.restype = ctypes.c_double
		fortran.krome_get_mu.argtypes = [array_1d_double]
		## extern char **krome_get_rnames();
		fortran.krome_get_mass.restype = array_1d_double
		fortran.krome_get_mass.argtypes = None
		fortran.krome_get_imass.restype = array_1d_double
		fortran.krome_get_imass.argtypes = None
		fortran.krome_get_hnuclei.restype = ctypes.c_double
		fortran.krome_get_hnuclei.argtypes = [array_1d_double]
		fortran.krome_get_charges.restype = array_1d_double
		fortran.krome_get_charges.argtypes = None
		## extern char **krome_get_names();
		##extern int krome_get_index(char *name);
		fortran.krome_get_rho.restype = ctypes.c_double
		fortran.krome_get_rho.argtypes = [array_1d_double]
		fortran.krome_scale_z.restype = None
		fortran.krome_scale_z.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_set_z.restype = None
		fortran.krome_set_z.argtypes = [ctypes.c_double]
		fortran.krome_set_clump.restype = None
		fortran.krome_set_clump.argtypes = [ctypes.c_double]
		fortran.krome_get_electrons.restype = ctypes.c_double
		fortran.krome_get_electrons.argtypes = [array_1d_double]
		fortran.krome_print_best_flux.restype = None
		fortran.krome_print_best_flux.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int]
		fortran.krome_print_best_flux_frac.restype = None
		fortran.krome_print_best_flux_frac.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double]
		fortran.krome_print_best_flux_spec.restype = None
		fortran.krome_print_best_flux_spec.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]
		fortran.krome_get_flux.restype = array_1d_double
		fortran.krome_get_flux.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_explore_flux.restype = None
		fortran.krome_explore_flux.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]
		fortran.krome_get_qeff.restype = array_1d_double
		fortran.krome_get_qeff.argtypes = None
		#IFKROME_useStars
		# Here is another example of a Fortran subroutine which accepts optional arguments:
		fortran.krome_stars_coe.restype = array_1d_double
		fortran.krome_stars_coe.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, array_1d_double, array_1d_double]
		fortran.krome_stars_energy.restype = array_1d_double
		fortran.krome_stars_energy.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_double, array_1d_double]
		#ENDIFKROME
		fortran.krome_dump_flux.restype = None
		fortran.krome_dump_flux.argtypes = [array_1d_double, ctypes.c_double, ctypes.c_int]
		fortran.krome_dump_rates.restype = None
		fortran.krome_dump_rates.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]
		fortran.krome_get_info.restype = None
		fortran.krome_get_info.argtypes = [array_1d_double, ctypes.c_double]
		fortran.krome_get_jacobian.restype = None
		fortran.krome_get_jacobian.argtypes = [ctypes.c_int, array_1d_double, ctypes.c_double, array_1d_double]

