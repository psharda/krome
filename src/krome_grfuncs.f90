!This module contains functions and subroutines
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.
module krome_grfuncs
contains

#KROME_header

  !**********************
  !get Tdust from tables, K
  function get_table_Tdust(n) result(Tdust)
    use krome_commons
    use krome_fit
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::ntot,Tdust,Tgas

    Tgas = n(idx_Tgas)

    !default, K
    Tdust = 1d0

    !total densitym, cm-3
    ntot = sum(n(1:nmols))

    !zero density returns default
    if(ntot==0d0) return

#IFKROME_dust_table_2D
    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
         dust_mult_Tgas, &
         log10(ntot), log10(Tgas))
#ENDIFKROME

#IFKROME_dust_table_3D
    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab3D(dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_AvVariable(:), &
         dust_tab_Tdust(:,:,:), dust_mult_ngas, &
         dust_mult_Tgas, dust_mult_AvVariable, &
         log10(ntot), log10(Tgas), dust_table_AvVariable_log)
#ENDIFKROME

  end function get_table_Tdust

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
         * pre_kvgas_sqrt * ims * sqrTgas &
         * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
         + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
       expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
         + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !******************
  function krate_2bodySi(n,idx1,idx2,Ea,Tdust) result(krate)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec),Ea,Tdust
    integer,intent(in)::idx1,idx2
    real*8::krate,amin,amax,pexp,d2g,rho0

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    krate = krate_2body(n(:),idx1,idx2,amin,amax,pexp,d2g,rho0,Ea,Tdust)

  end function krate_2bodySi

  !********************
  function krate_2body(n,idx1,idx2,amin,amax,pexp,d2g,rho0, &
       Ea,Tdust) result(krate)
    use krome_commons
    use krome_constants
    use krome_getphys
    implicit none
    integer,intent(in)::idx1,idx2
    real*8,intent(in)::n(nspec),amin,amax,pexp,d2g,rho0,Ea,Tdust
    real*8::rhog,p3,p4,ndns,krate,mred,fice,fbare,Preac
    real*8::iTd23,Ebare(nspec),Eice(nspec),mass(nspec)
    real*8,parameter::app2=(3d-8)**2 !cm^2 (Hocuk+2015)
    real*8,parameter::nu0=1d12 !1/s
    real*8,parameter::hbar=planck_erg/2d0/pi !erg*s
    real*8,parameter::ar=1d-8 !cm

    mass(:) = get_mass()

    !gas density, g/cm3
    rhog = sum(mass(1:nmols)*n(1:nmols))

    !exponentes
    p3 = pexp + 3d0
    p4 = pexp + 4d0

    !number of sites cm-3/mly
    ndns = rhog/(4d0/3d0*rho0*app2)*(amax**p3-amin**p3) &
         / (amax**p4-amin**p4) * p4 / p3

    !ice/bare fraction
    fbare = 1d0
#IFKROME_hasH2O
    fice = (n(idx_H2O_total)-n(idx_H2O))/ndns
    fbare = 1d0 - fice
#ENDIFKROME

    !reduced mass
    mred = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2))

    !tunneling probability
    Preac = exp(-2d0*ar/hbar*sqrt(2d0*mred*Ea*boltzmann_erg))

    !exponent
    iTd23 = 2d0/3d0/Tdust

    !get Ebind, K
    Ebare(:) = get_Ebind_bare()
#IFKROME_hasH2O
    Eice(:) = get_Ebind_ice()
#ENDIFKROME

    !compute rate
    krate = fbare*(exp(-Ebare(idx1)*iTd23)+exp(-Ebare(idx2)*iTd23))
#IFKROME_hasH2O
    krate = krate + fice*(exp(-Eice(idx1)*iTd23)+exp(-Eice(idx2)*iTd23))
#ENDIFKROME

    !rate in cm3/s
    krate = nu0*Preac/ndns*krate

  end function krate_2body

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
       dust_get_inv_phi(i) = 0d0
       dephi = (4d0 * nndust(i) * pi * asize2(i))
       if(dephi.le.0d0) cycle
       dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

#IFKROME_useChemisorption
  !***************************
  function dust_get_rateChem_PC(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_PC(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_PC(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_PC(idx+1)-dust_rateChem_PC(idx)) &
            + dust_rateChem_PC(idx)
    end do

  end function dust_get_rateChem_PC

  !***************************
  function dust_get_rateChem_CP(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_CP(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_CP(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_CP(idx+1)-dust_rateChem_CP(idx)) &
            + dust_rateChem_CP(idx)
    end do

  end function dust_get_rateChem_CP

  !***************************
  function dust_get_rateChem_CC(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_CC(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_CC(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_CC(idx+1)-dust_rateChem_CC(idx)) &
            + dust_rateChem_CC(idx)
    end do

  end function dust_get_rateChem_CC
#ENDIFKROME

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
       Td100 = Tdust(i) * 1d-2
       dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
            + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

 !*************************
  function dust_stick(Tgas,Tdust)
    implicit none
    real*8,intent(in)::Tgas,Tdust
    real*8::dust_stick
    real*8::Tg100,Td100

      Tg100 = Tgas * 1d-2
      Td100 = Tdust * 1d-2
      dust_stick = 1d0/(1d0 + 0.4d0*sqrt(Tg100+Td100) &
           + 0.2d0*Tg100 + 0.08d0*Tg100**2)

  end function dust_stick

  !****************************
  !sticking rate (1/s), assuming power-law dust distribution
  ! example rate is
  !  @format:idx,R,P,rate
  !  1,CO,CO_ice,krate_stick(n(:),idx_CO,1d-7,1d-5,-3.5,3d0,1d-2)
  ! n(:): internal status array (number densities, temeperature, etc...)
  ! idx : index of the sticking species, e.g. idx_CO
  ! Tdust: dust temperature (assume same for all bins), K
  ! amin: min grain size, cm
  ! amax: max grain size, cm
  ! pexp: power-law exponent, usually -3.5
  ! rho0: bulk material density, g/cm3, e.g. 3 g/cm3 for silicates
  ! d2g: dust to gass mass ratio, usually 0.01
  function krate_stick(n,idx,Tdust,amin,amax,pexp,rho0,d2g) result(k)
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    real*8,intent(in)::n(nspec),Tdust,amin,amax,pexp,rho0,d2g
    real*8::k,imass(nspec),p4,p3,mass(nspec),rhod
    integer,intent(in)::idx

    !get inverse mass squared
    imass(:) = get_imass_sqrt()
    !get masses
    mass(:) = get_mass()
    !derived exponents
    p3 = pexp + 3.
    p4 = pexp + 4.

    !total dust density, g/cm3
    rhod = sum(n(1:nmols)*mass(1:nmols))*d2g

    !compute rate (1/s) coefficient assuming normalization
    k = pre_kvgas_sqrt*sqrt(n(idx_Tgas)) * imass(idx) &
         * rhod / (4./3.*rho0) * p4 / p3 &
         * (amax**p3-amin**p3) / (amax**p4-amin**p4) &
         * dust_stick(n(idx_Tgas),Tdust)

  end function krate_stick

  !********************************
  !compact version of krate_stick
  function krate_stickSi(n,idx,Tdust) result(k)
    use krome_commons
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,amin,amax,d2g,rho0,pexp

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    k = krate_stick(n(:),idx,Tdust,amin,amax,pexp,rho0,d2g)

  end function krate_stickSi

  !***************************
  !evaporation rate, 1/s
  function krate_evaporation(n,idx,Tdust) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,Ebind(nspec),nu0

    nu0 = 1d12 !1/s
    Ebind(:) = get_EbindBare()

    k = nu0 * exp(-Ebind(idx)/Tdust)

  end function krate_evaporation

  !***************************
  !non-thermal evaporation rate (1/s) following Hollenbach 2009,
  ! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0809.1642
  !Gnot is the habing flux (1.78 is Draine)
  !Av is the visual extinction
  !crflux the ionization flux of cosmic rays, 1/s
  !yield is the efficiency of the photons to desorb the given molecule
  function krate_nonthermal_evaporation(idx, Gnot, Av, crflux, yield) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,parameter::crnot=1.3d-17
    real*8,parameter::Fnot=1d8 !desorbing photons flux, 1/s
    real*8,parameter::ap2=(3d-8)**2 !sites separation squared, cm2
    real*8,intent(in)::Gnot, Av, crflux, yield
    real*8::k,f70,kevap70(nspec)

    f70 = 3.16d-19*crflux/crnot
    kevap70(:) = get_kevap70()

    k = Gnot*Fnot*ap2*yield*exp(-1.8*Av)
    k = k + f70*kevap70(idx)

  end function krate_nonthermal_evaporation

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
       dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

#KROME_Ebareice

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

#KROME_Ebareice23

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

#KROME_Ebind_ice

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

#KROME_Ebind_bare

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

#KROME_parent_dust_bin

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table


  !**********************
  ! Cluster growth rate based on kinetic nucleation theory (KNT)
  ! Theory is explained in chapter 13 of Gail and Sedlmayr 2013
  ! (https://doi.org/10.1017/CBO9780511985607)
  function cluster_growth_rate(monomer_idx, cluster_size, temperature, stick) result(rate)
    ! k_N = v_thermal * cross_sectrion_N * stick_N
    ! with N the cluster size of the reactant
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    integer, parameter :: dp=kind(0.d0) ! double precision

    integer, intent(in) :: monomer_idx
    integer, intent(in) :: cluster_size
    real(dp), intent(in) :: temperature
    real(dp), intent(in), optional :: stick
    real(dp) :: rate

    real(dp) :: v_thermal
    real(dp) :: cross_section
    real(dp) :: stick_coefficient
    real(dp) :: monomer_radius
    real(dp) :: monomer_mass
    real(dp) :: mass(nspec)

    mass(:) = get_mass()

    if(monomer_idx == idx_TiO2) then
      ! Interatomic distance from Jeong et al 2000 DOI:10.1088/0953-4075/33/17/319
      monomer_radius = 1.78e-8_dp ! in cm
    else
      print *, "Monomer radius not yet defined"
    end if

    monomer_mass = mass(monomer_idx)

    v_thermal = sqrt(8._dp * boltzmann_erg * temperature &
              / (pi * monomer_mass))

    cross_section = pi * monomer_radius**2._dp * cluster_size**(2._dp/3._dp)

    ! Sticking coefficiet is set to one for simplicity
    if(present(stick)) then
      stick_coefficient = stick
    else
      stick_coefficient = 1._dp
    end if

    rate = v_thermal * cross_section * stick_coefficient

  end function cluster_growth_rate


  !**********************
  ! Cluster destruction rate based on kinetic nucleation theory (KNT)
  ! Theory is explained in chapter 13 of Gail and Sedlmayr 2013
  ! (https://doi.org/10.1017/CBO9780511985607)
  ! This reversed reaction is infered from detailed balance
  function cluster_destruction_rate(monomer_idx, n, cluster_size,&
     temperature, stick) result(rate)
    ! k_N = v_thermal * cross_section_(N-1) * stick_(N-1)
    ! * [n_1 * n_(N-1)/n_N]_equilibrium
    ! with N the cluster size of the reactant
    ! and [n_1 * n_(N-1)/n_N]_equilibrium are numbers densities in equilibrium
    ! k_N_destr = k_(N-1)_growth * [n_1 * n_(N-1)/n_N]_equi
    use krome_constants
    use krome_commons
    implicit none
    integer, parameter :: dp=kind(0.d0) ! double precision

    integer, intent(in) :: monomer_idx
    real(dp), intent(in) :: n(nspec)
    integer, intent(in) :: cluster_size
    real(dp), intent(in) :: temperature
    real(dp), intent(in), optional :: stick
    real(dp) :: rate

    real(dp) :: k_growth, ngas, pressure_scaled
    real(dp) :: gibbs_big, gibbs_small, gibbs_monomer
    real(dp) :: gibbs_part, gibbs_corr

    ! [n_(N-1)/n_N]_equi = (n_gas/n_1_equi) * exp( (dG_N - dG_(N-1) - dG_1) / RT )
    ! with dG_N "is the change in free enthalpy in the reaction of formation
    ! of 1 mol of clusters of size N from N mol of monomers."
    ! - Gail & Sedlmayr 2013, sec. 13.4.1
    ! See Clouet 2010 https://arxiv.org/abs/1001.4131v2
    ! pages 15+ for a more detailed derivation of
    ! the Gibbs free enegery of the system.
    ! This assumes the clusters to be dilute compared to the total gas
    ! which is resonable


    ! ngas = everything besides the clusters, but as clusters are assumed to be dilute
    ! their number density can be neglected compared to the total
    ! ngas = sum(n(1:nmols))
    if(monomer_idx == idx_TiO2 )then
      ngas = sum(n(1:nmols))-sum(n(idx_TiO2:idx_Ti10O20))
    else
      print *, "Clusters other than TiO2 are not yet defined"
    endif

    ! total gas pressure in units of 1 bar
    pressure_scaled = ngas * boltzmann_erg * temperature * 1.e-6_dp

    gibbs_big = gibbs_free_energy(monomer_idx, cluster_size, temperature) ! kJ*mol**(-1)
    gibbs_small = gibbs_free_energy(monomer_idx, cluster_size-1, temperature)! kJ*mol**(-1)
    gibbs_monomer = gibbs_free_energy(monomer_idx, 1, temperature)! kJ*mol**(-1)
    ! correction to the Gibbs free enegery under non-standard pressure of 1 bar.
    ! This only differs in the translational partition function.
    gibbs_corr = temperature * Rgas_kJ *log(pressure_scaled)

    ! gibss correction needs to be added to each gibss energy but _big and _small cancel
    gibbs_part = exp( (gibbs_big - gibbs_small - gibbs_monomer - gibbs_corr)&
                / ( Rgas_kJ * temperature ) )

    k_growth = cluster_growth_rate(monomer_idx, cluster_size-1, temperature)



    rate = k_growth * ngas * gibbs_part ! s^(-1)

  end function cluster_destruction_rate


  !**********************
  ! Change in free enthalpy in the reaction of formation
  ! of 1 mol of clusters of size N from N mol of monomers."
  ! - Gail & Sedlmayr 2013, sec. 13.4.1
  function gibbs_free_energy(monomer_idx, cluster_size, temperature) result(gibbs)
    use krome_constants
    use krome_commons
    implicit none
    integer, parameter :: dp=kind(0.d0) ! double precision

    integer, intent(in) :: monomer_idx
    integer, intent(in) :: cluster_size
    real(dp), intent(in) :: temperature
    real(dp) :: gibbs

    real(dp) :: Tinv, T, T2, T3
    real(dp) :: a, b, c, d, e

    T = temperature
    Tinv = T**(-1._dp)
    T2 = T*T
    T3 = T2*T

    if(monomer_idx == idx_TiO2) then
      ! Data taken from Lee, G et al. 2015 (10.1051/0004-6361/201424621)
      ! Lee, G. et al. 2018 fitted their own results (https://arxiv.org/pdf/1801.08482.pdf)
      ! dG = a*T**-1 + b + c*T + d*T**2 + e*T**3
      ! This fit is valid for 500 < T < 2000 K
      if(cluster_size == 1) then
        a = -1.63472903e3_dp
        b = -2.29197239e2_dp
        c = -3.60996766e-2_dp
        d = 1.60056318e-5_dp
        e =-2.02075337e-9_dp
      else if(cluster_size == 2) then
        a = -4.39367806e3_dp
        b = -9.77431160e2_dp
        c = 1.01656231e-1_dp
        d = 2.16685151e-5_dp
        e = -2.90960794e-9_dp
      else if(cluster_size == 3) then
        a = -7.27464297e3_dp
        b = -1.72789122e3_dp
        c = 2.40409836e-1_dp
        d = 2.74002833e-5_dp
        e = -3.81294573e-9_dp
      else if(cluster_size == 4) then
        a = -1.02808569e4_dp
        b = -2.51074121e3_dp
        c = 4.15061961e-1_dp
        d = 3.30076021e-5_dp
        e = -4.69138304e-9_dp
      else if(cluster_size == 5) then
        a = -1.37139638e4_dp
        b = -3.27506794e3_dp
        c = 5.73212328e-1_dp
        d = 4.12461166e-5_dp
        e = -6.14829810e-9_dp
      else if(cluster_size == 6) then
        a = -1.60124756e4_dp
        b = -4.13772573e3_dp
        c = 7.32672450e-1_dp
        d = 4.44131101e-5_dp
        e = -6.48290229e-9_dp
      else if(cluster_size == 7) then
        a = -1.89334054e4_dp
        b = -4.91964308e3_dp
        c = 8.93689186e-1_dp
        d = 4.99942488e-5_dp
        e = -7.35905348e-9_dp
      else if(cluster_size == 8) then
        a = -2.17672541e4_dp
        b = -5.72492348e3_dp
        c = 1.05703014e0_dp
        d = 5.57819924e-5_dp
        e = -8.27043313e-9_dp
      else if(cluster_size == 9) then
        a = -2.48377680e4_dp
        b = -6.51357184e3_dp
        c = 1.22288686e0_dp
        d = 6.10116309e-5_dp
        e = -0.08225913e-9_dp
      else if(cluster_size == 10) then
        a = -2.76078426e4_dp
        b = -7.34516329e3_dp
        c = 1.37500651e0_dp
        d = 6.70631142e-5_dp
        e = -1.00410219e-8_dp
      else
        print *, "There is no thermochemical data on &
         TiO2 clusters larger than 10."
       end if
       gibbs = a*Tinv + b + c*T + d*T2 + e*T3 ! kJ*mol**(-1)

     else
       print *, "There is no thermochemical data on &
       clusters ofther than TiO2."
    end if

  end function gibbs_free_energy

end module krome_grfuncs
