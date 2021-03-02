# **************************
# Draine1978 flux eV/cm2/sr
def Jdraine(energy):
    planck_eV = 4.135667662e-15  # eV*s
    if energy > 13.6 or energy < 5e0:
        return 0e0
    return (1.658e6 * energy - 2.152e5 * energy ** 2 + 6.919e3 * energy ** 3) * energy * planck_eV


# **************************
# integral of Draine, J(E)/E
def intJdraine():
    planck_eV = 4.135667662e-15  # eV*s
    energyMin = 6.1992  # ev, 200 nm
    energyMax = 13.6

    def f(energy):
        return (1.658e6 * energy ** 2 / 2e0 - 2.152e5 * energy ** 3 / 3e0 + 6.919e3 * energy ** 4 / 4e0) * planck_eV

    return f(energyMax) - f(energyMin)


# **************************
# black body radiance @ Tbb, eV/cm2/sr
def JBB(energy, Tbb):
    from math import exp
    planck_eV = 4.135667662e-15  # eV*s
    boltzmann_eV = 8.6173303e-5  # eV/K
    clight = 2.99792458e10  # cm/s

    # exponential arg
    xexp = energy / boltzmann_eV / Tbb

    # limit exp overflow
    if xexp > 3e2 or energy < 1e-10: return 0e0
    return 2e0 * energy ** 3 / planck_eV ** 2 / clight ** 2 \
           / (exp(xexp) - 1e0)


# ************************
# integral of black body radiance @ Tbb, J(E)/E
def JBBe(energy, Tbb):
    return JBB(energy, Tbb) / energy


# **************************
# integrate black body J(E)/E
def intJBB(Tbb):
    from scipy.integrate import quad
    return quad(JBBe, 6.1992, 13.6, args=(Tbb), epsabs=1e-20)[0]
