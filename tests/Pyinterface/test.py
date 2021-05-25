# -*- coding: utf-8 -*-
import numpy as np
from pykrome import PyKROME
import ctypes

if __name__ == "__main__":

    pyk = PyKROME()
    pyk.lib.krome_init()

    spy = pyk.krome_seconds_per_year

    # ----==== Isotope decay off ====---- #

    x = np.ones(pyk.krome_nmols) * 1e-20
    ntot = 1.0e4 # cm**-3
    x[pyk.krome_idx_Hj] = ntot

    # set iron abundance
    pyk.lib.krome_scale_z(x, 0.0)
    x[pyk.krome_idx_FEj] = x[pyk.krome_idx_FE]
    x[pyk.krome_idx_FE] = 0.0
    x[pyk.krome_idx_E] = pyk.lib.krome_get_electrons(x)

    # set 60Fe abundance
    pyk.lib.krome_set_user_tauh(1.5e6 * spy)
    x_Fe60 = 0.0
    print('Fe-60 abundance =', x_Fe60)
    x[pyk.krome_idx_60FE] = x_Fe60 * x[pyk.krome_idx_FEj]
    # set rate for ionisation by isotope decay
    fe60_xi = 0.0
    print('Isotope decay off.')
    pyk.lib.krome_set_user_xi(fe60_xi)
    # set heating from isotope decay
    pyk.lib.krome_set_user_wergs(36.0 * pyk.krome_eV_to_erg)

    # in order to have the updated `Tgas` that KROME returns tracked, `Tgas`
    # needs to be a ctype.
    Tgas = ctypes.c_double(1.0e3)

    t = 0.0
    dt = 1.0e-2 * spy
    nstep = 0
    output = []
    while t <= 1.0e8 * spy:
        dt = dt * 1.1
        t += dt
        if np.mod(nstep,10) == 0:
          print("nstep = {0:4d}".format(nstep))

        # call KROME
        pyk.lib.krome(x, ctypes.byref(Tgas), ctypes.byref(ctypes.c_double(dt)))
        nstep += 1

        output.append(np.concatenate((np.array([t/spy, Tgas.value]), x/ntot)))

        pyk.lib.krome_popcool_dump(t/spy, 71)

    # write output
    output = np.array(output)
    np.savetxt('idoff.py.dat',output,fmt='%15.8E',delimiter='')

    print("Finished. Number of steps = {}".format(nstep))

    # ----==== Isotope decay on ====---- #

    x = np.ones(pyk.krome_nmols) * 1e-20
    ntot = 1.0e4 # cm**-3
    x[pyk.krome_idx_Hj] = ntot

    # set iron abundance
    pyk.lib.krome_scale_z(x, 0.0)
    x[pyk.krome_idx_FEj] = x[pyk.krome_idx_FE]
    x[pyk.krome_idx_FE] = 0.0
    x[pyk.krome_idx_E] = pyk.lib.krome_get_electrons(x)

    # set 60Fe abundance
    pyk.lib.krome_set_user_tauh(1.5e6 * spy)
    x_Fe60 = 1.0e-6
    print('Fe-60 abundance =', x_Fe60)
    x[pyk.krome_idx_60FE] = x_Fe60 * x[pyk.krome_idx_FEj]
    # set rate for ionisation by isotope decay
    fe60_xi = 1.0e-10
    print('Isotope decay on.')
    pyk.lib.krome_set_user_xi(fe60_xi)
    # set heating from isotope decay
    pyk.lib.krome_set_user_wergs(36.0 * pyk.krome_eV_to_erg)

    # in order to have the updated `Tgas` that KROME returns tracked, `Tgas`
    # needs to be a ctype.
    Tgas = ctypes.c_double(1.0e3)

    t = 0.0
    dt = 1.0e-2 * spy
    nstep = 0
    output = []
    while t <= 1.0e8 * spy:
        dt = dt * 1.1
        t += dt
        if np.mod(nstep, 10) == 0:
          print("nstep = {0:4d}".format(nstep))

        # call KROME
        pyk.lib.krome(x, ctypes.byref(Tgas), ctypes.byref(ctypes.c_double(dt)))
        nstep += 1

        output.append(np.concatenate((np.array([t/spy, Tgas.value]), x/ntot)))

        pyk.lib.krome_popcool_dump(t/spy, 72)

    # write output
    output = np.array(output)
    np.savetxt('idon.py.dat',output,fmt='%15.8E',delimiter='')

    print("Finished. Number of steps = {}".format(nstep))
