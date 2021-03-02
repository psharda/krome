############## rate functions #############
# This file contains special rate funtions
# that are internally defined in KROME
# but not in the network file

# List of all special functions
functionList = ["cluster_growth_rate",
                "cluster_destruction_rate",
                "steady_state_nucleation_rate",
                ]


def cluster_growth_rate(monomer, cluster_size, temperature, stick=1.0):
    # import constants
    # TODO: make file with constants
    from math import sqrt
    pi = 3.14159265359e0
    boltzmann_erg = 1.380648e-16  # erg / K

    # if monomer_idx == "idx_TiO2":
    #     # Interatomic distance from Jeong et al 2000 DOI:10.1088/0953-4075/33/17/319
    #     monomer_radius = 1.78e-8  # in cm
    #     monomer_name = "TIO2"
    # else:
    #     print "Monomer radius not yet defined"

    monomer_radius = monomer.radius
    # monomer_radius = 1.78e-8  # in cm
    inv_monomer_mass = 1. / monomer.mass
    inv_cluster_mass = 1. / cluster_size * inv_monomer_mass
    inv_reduced_mass = inv_monomer_mass + inv_cluster_mass

    v_thermal = sqrt(8.0 * boltzmann_erg * temperature * inv_reduced_mass
                     / pi)
    cluster_radius = monomer_radius * cluster_size ** (1. / 3.)
    cross_section = pi * (monomer_radius + cluster_radius) ** 2.

    rate = v_thermal * cross_section * stick

    return rate


# **********************
# Cluster destruction rate based on kinetic nucleation theory (KNT)
# Theory is explained in chapter 13 of Gail and Sedlmayr 2013
# (https://doi.org/10.1017/CBO9780511985607)
# This reversed reaction is infered from detailed balance
# NOTE: n is currently fixed (for testing)###
def cluster_destruction_rate(monomer, cluster_size,
                             temperature, stick=1.):
    # k_N = v_thermal * cross_section_(N-1) * stick_(N-1)
    # * [n_1 * n_(N-1)/n_N]_equilibrium
    # with N the cluster size of the reactant
    # and [n_1 * n_(N-1)/n_N]_equilibrium are numbers densities in equilibrium
    # k_N_destr = k_(N-1)_growth * [n_1 * n_(N-1)/n_N]_equi
    from math import log, exp
    boltzmann_erg = 1.380648e-16  # erg / K
    Rgas_kJ = 8.3144621e-3  # kJ/K/mol
    # [n_(N-1)/n_N]_equi = (n_gas/n_1_equi) * exp( (dG_N - dG_(N-1) - dG_1) / RT )
    # with dG_N "is the change in free enthalpy in the reaction of formation
    # of 1 mol of clusters of size N from N mol of monomers."
    # - Gail & Sedlmayr 2013, sec. 13.4.1
    # See Clouet 2010 https://arxiv.org/abs/1001.4131v2
    # pages 15+ for a more detailed derivation of
    # the Gibbs free enegery of the system.
    # This assumes the clusters to be dilute compared to the total gas
    # which is resonable

    ### uncorrected rate
    # ngas = everything besides the clusters, but as clusters are assumed to be dilute
    # their number density can be neglected compared to the total
    # ngas = sum(n(1:nmols))
    # total gas number density
    # ngas = n # cm^(-3)
    ###

    gibbs_big = gibbs_free_energy(monomer, cluster_size, temperature)  # kJ*mol**(-1)
    gibbs_small = gibbs_free_energy(monomer, cluster_size - 1, temperature)  # kJ*mol**(-1)
    gibbs_monomer = gibbs_free_energy(monomer, 1, temperature)  # kJ*mol**(-1)
    # correction to the Gibbs free enegery under non-standard pressure of 1 bar.
    # This only differs in the translational partition function.
    # # total gas pressure in units of 1 bar
    # pressure_scaled = ngas * boltzmann_erg * temperature * 1.e-6
    # gibbs_corr = temperature * Rgas_kJ * log(pressure_scaled)
    # gibss correction needs to be added to each gibss energy but _big and _small cancel
    # The gibbs_corr factor ultimately cancels out ngas and reduced to:
    non_standard_correction = (1.e-6 * boltzmann_erg * temperature) ** (-1)

    gibbs_part = (exp((gibbs_big - gibbs_small - gibbs_monomer)
                      / (Rgas_kJ * temperature)))

    k_growth = cluster_growth_rate(monomer, cluster_size - 1, temperature)

    # uncorrected rate
    # rate = k_growth * ngas * gibbs_part
    # corrected rate
    rate = k_growth * gibbs_part * non_standard_correction  # s^(-1)

    return rate


# nucleation rate based on the assumption of a steady state (number densities
# of the cluster don't change anymore). More details on how this is derived:
# see chapter 13.3.5 of Gail and Sedlmayr 2013
# (https://doi.org/10.1017/CBO9780511985607)
# NOTE: this is not useful when evolving a kinetic network
# which is inhernetly time dependent.
# However, it is frequently used in literature, and therefore useful as
# a way of comparing (mainly qualitatively/roughly)
def steady_state_nucleation_rate(monomer, ngas, max_cluster_size, temperature):
    # ngas can be changed with a different input parameter
    # e.g. total or partial pressure
    boltzmann_erg = 1.380648e-16
    p_mass = 1.67262158e-24
    monomer_frac = 1e-10
    # Frequently used units (Who even thinks in pressure?! o.O)
    # convert pressure (in bar) to number density (in cm^(-3))
    # ngas = pressure/(boltzmann_erg * temperature * 1.e-6)
    # convert 1 bar tonumber density (in cm^(-3))
    # ngas = 1. /(boltzmann_erg * temperature * 1.e-6)
    # convert mass density (in g cm^(-3)) to number density (in cm^(-3))
    # mu = 2.35 (70/30 H2/He in mass density)
    # ngas = rhogas / (mu * p_mass)

    n_monomer = ngas * monomer_frac
    # or get directly/conversion as input
    # n_monomer = pressure/(boltzmann_erg * temperature)

    # See eq. 13.17 for definition of rate with last term << 1 (in s^-1 cm^-3 )
    # We prefer to work with the growth and destruction rate as described
    # in above functions (see Eq. II in Bromley+2016 DOI:10.1039/c6cp03629e)
    part1 = cluster_growth_rate(monomer, 1, temperature) * n_monomer ** 2

    sumpart = 0.
    for j in range(2, max_cluster_size + 1):
        kmin = 1.
        kplus = 1.
        for i in range(2, j + 1):
            kmin *= cluster_destruction_rate(monomer, i, temperature)
            kplus *= n_monomer * cluster_growth_rate(monomer, i - 1, temperature)

        sumpart += (kmin / kplus)

    totpart = 1 + sumpart

    rate = part1 / totpart

    return rate


# **********************
# Change in free enthalpy in the reaction of formation
# of 1 mol of clusters of size N from N mol of monomers."
# - Gail & Sedlmayr 2013, sec. 13.4.1
def gibbs_free_energy(monomer, cluster_size, temperature):
    T = temperature
    Tinv = T ** (-1.)
    T2 = T * T
    T3 = T2 * T

    if monomer.name == "TiO2":
        # Data taken from Lee, G et al. 2015 (10.1051/0004-6361/201424621)
        # Lee, G. et al. 2018 fitted their own results (https://arxiv.org/pdf/1801.08482.pdf)
        # dG = a*T**-1 + b + c*T + d*T**2 + e*T**3
        # This fit is valid for 500 < T < 2000 K
        if cluster_size == 1:
            a = -1.63472903e3
            b = -2.29197239e2
            c = -3.60996766e-2
            d = 1.60056318e-5
            e = -2.02075337e-9
        elif cluster_size == 2:
            a = -4.39367806e3
            b = -9.77431160e2
            c = 1.01656231e-1
            d = 2.16685151e-5
            e = -2.90960794e-9
        elif cluster_size == 3:
            a = -7.27464297e3
            b = -1.72789122e3
            c = 2.40409836e-1
            d = 2.74002833e-5
            e = -3.81294573e-9
        elif cluster_size == 4:
            a = -1.02808569e4
            b = -2.51074121e3
            c = 4.15061961e-1
            d = 3.30076021e-5
            e = -4.69138304e-9
        elif cluster_size == 5:
            a = -1.37139638e4
            b = -3.27506794e3
            c = 5.73212328e-1
            d = 4.12461166e-5
            e = -6.14829810e-9
        elif cluster_size == 6:
            a = -1.60124756e4
            b = -4.13772573e3
            c = 7.32672450e-1
            d = 4.44131101e-5
            e = -6.48290229e-9
        elif cluster_size == 7:
            a = -1.89334054e4
            b = -4.91964308e3
            c = 8.93689186e-1
            d = 4.99942488e-5
            e = -7.35905348e-9
        elif cluster_size == 8:
            a = -2.17672541e4
            b = -5.72492348e3
            c = 1.05703014e0
            d = 5.57819924e-5
            e = -8.27043313e-9
        elif cluster_size == 9:
            a = -2.48377680e4
            b = -6.51357184e3
            c = 1.22288686e0
            d = 6.10116309e-5
            e = -0.08225913e-9
        elif cluster_size == 10:
            a = -2.76078426e4
            b = -7.34516329e3
            c = 1.37500651e0
            d = 6.70631142e-5
            e = -1.00410219e-8
        else:
            print("There is no thermochemical data on "
                  "TiO2 clusters larger than 10.")

        gibbs = a * Tinv + b + c * T + d * T2 + e * T3  # kJ*mol**(-1)

    else:
        print("There is no thermochemical data on "
              "clusters ofther than TiO2.")

    return gibbs
