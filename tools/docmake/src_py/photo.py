def Jdraine(energy):
	planck_eV = 4.135667662e-15 #eV*s
	if(energy>13.6 or energy<5e0): return 0e0
	return (1.658e6*energy - 2.152e5*energy**2 + 6.919e3*energy**3) * energy * planck_eV
