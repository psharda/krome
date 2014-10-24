module krome_subs
contains

#KROME_header

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
#KROME_shortcut_variables
    real*8::small,nmax
    integer::i
#KROME_initcoevars
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = #KROME_small

#KROME_Tshortcuts

#KROME_coevars

    k(:) = small !inizialize coefficients

#KROME_krates

    coe(:) = k(:) !set coefficients to return variable


    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if

  end function coe

#KROME_metallicity_functions

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
       planckBB = 2d0*x**3/planck_eV**2/clight**2 &
            / (exp(xexp)-1d0)
    end if

  end function planckBB

  !****************************
  !tanh smoothing function that 
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
         + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that 
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
         + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0, 
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

#IFKROME_useShieldingDB96
    !compute shielding from Draine+Bertoldi 1996
    krome_fshield = calc_H2shieldDB96(n(:), Tgas)
#ENDIFKROME

#IFKROME_useShieldingWG11
    !compute shielding from Wolcott+Greene 2011
    krome_fshield =  calc_H2shieldWG11(n(:), Tgas)
#ENDIFKROME

  end function krome_fshield

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)
#KROME_conserve

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve


  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c  
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
         * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
         * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
         + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93


  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = n(idx_Tgas)

#KROME_num2col_method

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = n(idx_Tgas)

#KROME_col2num_method

  end function col2num

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
       fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
       fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2 
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
         0.035d0/sqrt(1d0+x) * &
         exp(-8.5d-4*sqrt(1+x))

  end function fselfH2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
       j1 = j + 1
       zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
       zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * zo**a*exp(-2d0*ibTgas)) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
       z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom 
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

#KROME_gamma

    gamma_index = krome_gamma
  end function gamma_index

  !*****************************
  !get the mean molecular weight in grams
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = sum(n(1:nmols)*m(1:nmols)) &
         / sum(n(1:nmols)) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight in grams
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / sum(n(1:nmols)) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

#KROME_masses

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

#KROME_imasses_sqrt

  end function get_imass_sqrt


  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

#KROME_imasses

  end function get_imass

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

#KROME_names

  end function get_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

#KROME_sum_H_nuclei
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

#KROME_zatoms

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

#KROME_qeff

  end function get_qeff

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
         /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
         /p_mass/gravity/mu)

  end function get_jeans_length_rho

#IFKROME_useShieldingDB96
  !************************
  !calculate the self-shielding factor, following Draine&Bertoldi 1996 
  !NOTE: this function is suited for collapse. Use with caution!
  function calc_H2shieldDB96(n,Tgas)
    use krome_commons
    real*8::n(nspec),Tgas,calc_H2shieldDB96,N_H2, nH2

    !check on H2 abundances to avoid
    ! weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))

    N_H2 = nH2*get_jeans_length(n(:),Tgas)*0.5d0  !column density (cm-2)
    calc_H2shieldDB96 = min(1.d0, (N_H2*1.d-14)**(-0.75d0))

  end function calc_H2shieldDB96
#ENDIFKROME

#IFKROME_useShieldingWG11
  !************************
  !calculate the self-shielding factor, following Wolcott&Greene 2011
  !NOTE: this function is suited for collapse. Use with caution!
  function calc_H2shieldWG11(n,Tgas)
    use krome_commons
    use krome_constants
    real*8::n(nspec),Tgas,calc_H2shieldWG11,N_H2,nH2
    real*8::xN_H2,b5,H_mass

    !check on H2 abundances to avoid weird numerical artifacts
    nH2 = max(1d-40, n(idx_H2))
    N_H2 = nH2*get_jeans_length(n(:) ,Tgas)*0.5d0  !column density (cm-2)
    xN_H2 = N_H2*2d-15 !normalized column density (#), 2d-15=1/5d14
    H_mass = p_mass+e_mass !H mass in g

    !doppler broadening parameter b divided by 1d5 cm/s (#)
    b5 = ((boltzmann_erg*Tgas/H_mass)**0.5d0)*1.d-5 
    calc_H2shieldWG11 = 0.965d0/(1.d0+xN_H2/b5)**1.1d0 &
         + (0.035d0/(1.d0+xN_H2)**0.5d0) &
         * exp(-8.5d-4*(1.d0+xN_H2)**0.5d0)

  end function calc_H2shieldWG11
#ENDIFKROME

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
         * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
         * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n, Tgas)
    use krome_commons
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
         69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
         5.313317d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
         -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
         0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
         21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
         2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
         -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
         0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
       logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
            k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
            k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
       logk_h2 = k_CIDm(7,i)*invT
       logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
            k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
       logk_l2 = k_CIDm(13,i)*invT      
       logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
            + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
       logn_c2 = k_CIDm(18,i) + logn_c1
       p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas*1.850d-3) &
            + k_CIDm(21,i)*exp(-Tgas*4.40d-2)
       n_c1 = 1d1**(-logn_c1)
       n_c2 = 1d1**(-logn_c2)
       logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H*n_c1)**p) &
            + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H*n_c2)**p)
       k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID 

  end function dissH2_Martin96

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nmol,nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nmol,nndust,ims,stick,adust2,sqrTgas

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

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
       a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
       exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

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


  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
       !when found store and break loop
       if(trim(names(i))== trim(name)) then
          get_index = i !store index
          exit
       end if
    end do

    !error if species not found
    if(get_index<0) then
       print *,"ERROR: can't find the index of ",name
       stop
    end if

  end function get_index

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

#KROME_electrons_balance
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

#KROME_charges

  end function get_charges

  !************************
  !get species charges
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

#KROME_reaction_names

  end function get_rnames

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  function revKc(Tgas,ridx,pidx)
    implicit none
    real*8::revKc,Tgas
    integer::ridx(:),pidx(:),i

    revKc = 0.d0

    do i=1,size(pidx)
       revKc = revKc + revHS(Tgas,pidx(i))
    end do

    do i=1,size(ridx)
       revKc = revKc - revHS(Tgas,ridx(i))
    end do

  end function revKc

  !*****************************
  !compute H-S for species with index idx 
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
#KROME_var_reverse
    integer::idx

    p1(:,:) = 0.d0
    p2(:,:) = 0.d0
    Tlim(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)

#KROME_kc_reverse

    if(Tlim(idx,2)==0.d0) then
       revHS = 0.d0
       return
    end if

    !select set of NASA polynomials using temperature
    if(Tlim(idx,1).le.Tgas .and. Tgas.le.Tlim(idx,2)) p(:) = p1(idx,:)
    if(Tlim(idx,2)<Tgas .and. Tgas.le.Tlim(idx,3)) p(:) = p2(idx,:)

    !compute NASA polynomials for enthalpy and enthropy
    H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
         p(5)*Tgas4*0.2d0 + p(6)*invT
    S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
         p(5)*Tgas4*0.25d0 + p(7)

    revHS = H - S

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes 
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
       print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes 
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
       print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
       stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
       if(flux(idx(i))<flux(idx(1))*frac) exit
       print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
       found = .false.
#KROME_arr_reactprod
       maxflux = max(maxflux,flux(i))
       if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
       print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
            flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local 

    n = size(f)
    !init indexes
    do i=1,n
       idx_sort(i) = i
    end do

    !loop to sort
    do
       found = .false. !swapped something flag
       do i=2,n
          !> for descending, < for ascending
          if(f(i)>f(i-1)) then
             found = .true.
             !swap real value
             ftmp = f(i)
             f(i) = f(i-1)
             f(i-1) = ftmp
             !swap index
             itmp = idx_sort(i)
             idx_sort(i) = idx_sort(i-1)
             idx_sort(i-1) = itmp
          end if
       end do
       !if nothing swapped exit
       if(.not.found) exit
    end do


  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
#KROME_rvars
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
#KROME_arrs
#KROME_arr_flux
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

#KROME_implicit_arrays

  end subroutine load_arrays

  !********************************************
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    character(len=*)::filename
    character(len=60)::row_string
    real*8::x(:),y(:),z(:,:),rout(3),xmul,ymul
    integer::i,j,ios

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
       print *,"ERROR: in init_anytab2D x size differs from z"
       stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
       print *,"ERROR: in init_anytab2D y size differs from z"
       stop
    end if

    print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(51,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
       stop
    end if

    !skip the comments and the first line with the sizes of the data
    ! which are already known from the pre-processing
    do
       read(51,'(a)') row_string
       if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
       print *,"ERROR: file "//filename//" should"
       print *," contain the number of rows and "
       print *," columns in the format"
       print *,"  RR, CC"
       print *,row_string
       stop
    end if

    !loop to read file
    do i=1,size(x)
       do j=1,size(y)
          read(51,*,iostat=ios) rout(:)
          y(j) = rout(2)
          z(i,j) = rout(3)
       end do
       x(i) = rout(1)
       read(51,*,iostat=ios) !skip blanks
       if(ios.ne.0) exit
    end do
    close(51)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !******************************
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j
    real*8::x(:),y(:),z(:,:),xmul,ymul,xx,yy,zz
    character(len=*)::fname

    open(91,file=fname//".fit",status="replace")
    open(92,file=fname//".org",status="replace")
    do i=1,size(x)
       do j=1,size(y)
          xx = x(i)
          yy = y(j)
          zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
          write(91,*) xx,yy,zz
          write(92,*) x(i),y(j),z(i,j)
       end do
       write(91,*)
       write(92,*)
    end do
    print *,"original file wrote in",fname//".org"
    print *,"fit test file wrote in",fname//".fit"

  end subroutine test_anytab2D

  !******************************
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  function fit_anytab1D(x,z,xmul,xx)
    real*8::fit_anytab1D,x(:),z(:),xmul,xx,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
       b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
       c(1) = 0d0
       d(1) = 0d0
       b(2) = b(1)
       c(2) = 0d0
       d(2) = 0d0
       return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
       d(i) = x(i+1) - x(i)
       b(i) = 2d0*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions 
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination 
    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
       ispline = y(1)
       return
    end if

    if(u>=x(n)) then
       ispline = y(n)
       return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
       k = (i+j)/2
       if(u<x(k)) then
          j=k
       else
          i=k
       end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_subs
