# -*- coding: utf-8 -*-
import numpy as np
import pykrome as pyk

if __name__ == "__main__":

    pyk.krome_init()

    spy = pyk.krome_seconds_per_year

    x = np.ones(pyk.krome_nmols) * 1e-20
    ntot = 1.0e4 # cm**-3
    x[pyk.krome_idx_Hj] = ntot

    # set iron abundance
    pyk.krome_scale_z(x, 0.0)
    x[pyk.krome_idx_Fej] = x[pyk.krome_idx_Fe]
    x[pyk.krome_idx_Fe] = 0.0
    x[pyk.krome_idx_E] = pyk.krome_get_electrons(x)

    # set 60Fe abundance
    pyk.krome_set_user_tauh(1.5e6 * spy)
    x[pyk.krome_idx_60Fe] = 1.0e-6 * x[pyk.krome_idx_Fej]
    # set rate for ionisation by isotope decay
    pyk.krome_set_user_xi(1.0e-10)
    # set heating from isotope decay
    pyk.krome_set_user_wergs(36.0 * pyk.krome_ev_to_erg)

    Tgas = 1.0e3

    t = 0.0
    dt = 1.0e-2 * spy
    nstep = 0
    output = []
    while t <= 1.0e8 * spy:
        dt = dt * 1.1
        t += dt
        print nstep, Tgas

        # call KROME
        Tgas = pyk.krome(x, Tgas, dt)
        nstep += 1

        output.append(np.concatenate((np.array([t/spy, Tgas]), x/ntot)))

        pyk.krome_popcool_dump(t/spy, 70)

    # write output
    output = np.array(output)
    np.savetxt('python.65',output,fmt='%17.8E',delimiter='')

    print "Finished. Number of steps = {}".format(nstep)
