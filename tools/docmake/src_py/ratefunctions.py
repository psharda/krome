############## rate functions #############
# This file contains special rate funtions
# that are internally defined in KROME
# but not in the network file

# List of all special functions
functionList = ["cluster_growth_rate",
                "cluster_destruction_rate",
                ]

def cluster_growth_rate(monomer, cluster_size, temperature, stick=1.0):
    # import constants
    from math import sqrt
    pi = 3.14159265359e0
    boltzmann_erg = 1.380648e-16 #erg / K

    # if monomer_idx == "idx_TiO2":
    #     # Interatomic distance from Jeong et al 2000 DOI:10.1088/0953-4075/33/17/319
    #     monomer_radius = 1.78e-8  # in cm
    #     monomer_name = "TIO2"
    # else:
    #     print "Monomer radius not yet defined"

    # monomer_radius = monomer.radius
    monomer_radius = 1.78e-8  # in cm
    monomer_mass = monomer.mass

    v_thermal = sqrt(8.0 * boltzmann_erg * temperature
                        / (pi * monomer_mass))
    cross_section = pi * monomer_radius**2 * cluster_size**(2./3.)

    rate = v_thermal * cross_section * stick

    return rate
