module krome_subs
contains

#KROME_header

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
#KROME_shortcut_variables
    real*8::small,nmax
    integer::i
#KROME_initcoevars
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

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


  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*50::fname,line
    integer::ios,i,nunit

    fname = "reactions_verbatim.dat"

#KROME_no_verbatim_file

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
       print *,"ERROR: "//trim(fname)//" file not present!"
       stop
    end if

    !load reactions from file
    do i=1,nrea
       read(nunit,'(a)',iostat=ios) line
       if(ios/=0) then
          print *,"ERROR: problem reading "//trim(fname)
          stop
       end if
       reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

#IFKROME_has_electrons
  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
         /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
         *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
         /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
         *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
         /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
         *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
         /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
         *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

#ENDIFKROME

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

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

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
#KROME_conserve_matrix
    B(:) = ref(:)

#IFKROME_useLAPACK
    call mydgesv(natoms,A(:,:),B(:), "conserveLin_x")
#ENDIFKROME

#KROME_conserve_fscale
#IFKROME_has_electrons
    !charge conservation
    x(idx_E) = m(idx_E)*(#KROME_conserveLin_electrons)
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
       print *,"ERROR in conserveLin, electrons < 0"
       stop
    end if
#ENDIFKROME

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

#KROME_conserveLin_ref

  end function conserveLinGetRef_x

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

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

#KROME_load_parts

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: partition function not found"
       print *," in file "//fname
       stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
       read(33,*,iostat=ios) rout(:)
       if(ios<0) exit
       if(ios.ne.0) cycle
       icount = icount + 1
       min_part = min(min_part,rout(1))
       array_tmp(icount) = rout(2)
       dT_part = rout(1) - Told
       Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part


  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
       kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
         + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

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

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
         * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

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
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
         69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
         5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
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
       p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
            + k_CIDm(21,i)*exp(-Tgas/4.40d2)
       n_c1 = 1d1**(logn_c1)
       n_c2 = 1d1**(logn_c2)
       logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
            + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
       k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96
#IFKROME_use_cluster_growth

  !**********************
  ! Cluster growth rate based on kinetic nucleation theory (KNT)
  ! Theory is explained in chapter 13 of Gail and Sedlmayr 2013
  ! (https://doi.org/10.1017/CBO9780511985607)
  function cluster_growth_rate(monomer_idx, cluster_size, temperature, stick) result(rate)
    ! k_N = v_thermal * cross_section_N * stick_N
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
    real(dp) :: cluster_radius
    real(dp) :: inverse_monomer_mass
    real(dp) :: inverse_cluster_mass
    real(dp) :: inverse_reduced_mass
    real(dp) :: inverse_mass(nspec)

    inverse_mass(:) = get_imass()
    
    !TODO: Don't hard code this here
    if(monomer_idx == idx_TiO2) then
      ! Interatomic distance from Jeong et al 2000 DOI:10.1088/0953-4075/33/17/319
      monomer_radius = 1.62e-8_dp ! in cm
    else if(monomer_idx == idx_Al2O3) then
      ! Interatomic distance O-Al-O (linear geometry) from Archibong et al 1999
      ! doi: 10.1021/jp983695n
      ! NOTE: temporary value
      monomer_radius = 3.304e-8_dp ! in cm
    else if(monomer_idx == idx_MgO) then
      ! Half a bond length form Farrow et al 2014 doi:10.1039/C4CP01825G
      monomer_radius = 0.865e-8_dp ! in cm
    else if(monomer_idx == idx_SiO) then
      ! Bond length from Bromley et al 2016 doi:10.1039/c6cp03629e
      monomer_radius = 0.75765e-8_dp ! in cm
    else
      print *, "Monomer radius not yet defined"
    end if

    inverse_monomer_mass = inverse_mass(monomer_idx)
    inverse_cluster_mass = 1._dp/cluster_size * inverse_monomer_mass
    inverse_reduced_mass = inverse_monomer_mass + inverse_cluster_mass

    v_thermal = sqrt(8._dp * boltzmann_erg * temperature &
              * inverse_reduced_mass / pi )

    ! Assuming cluster volume is proportional to monomer volume
    ! V_N = N * V_1, and both are considered as a hypothetical sphere
    cluster_radius = monomer_radius * cluster_size**(1._dp/3._dp)

    ! Geometrical cross section
    cross_section = pi * (monomer_radius + cluster_radius)**2._dp

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
  function cluster_destruction_rate(monomer_idx, cluster_size,&
     temperature, stick) result(rate)
    ! k_N = v_thermal * cross_section_(N-1) * stick_(N-1)
    ! * [n_1 * n_(N-1)/n_N]_equilibrium
    ! with N the cluster size of the reactant
    ! and [n_1 * n_(N-1)/n_N]_equilibrium are numbers densities in equilibrium
    ! k_N_destr = k_(N-1)_growth * [n_1 * n_(N-1)/n_N]_equi
    ! NOTE: this entire rate is equivalent to calling revKc() on
    ! cluster_growth_rate() given that the Gibbs free eneries (polynomials)
    ! of the clusters are stored in the thermochemical data of the species
    ! or provided as tables cfr. data/database/janaf/
    use krome_constants
    use krome_commons
    implicit none
    integer, parameter :: dp=kind(0.d0) ! double precision

    integer, intent(in) :: monomer_idx
    integer, intent(in) :: cluster_size
    real(dp), intent(in) :: temperature
    real(dp), intent(in), optional :: stick
    real(dp) :: rate

    real(dp) :: k_growth
    real(dp) :: gibbs_big, gibbs_small, gibbs_monomer
    real(dp) :: gibbs_part, non_standard_correction

    ! [n_(N-1)/n_N]_equi = (n_gas/n_1_equi) * exp( (dG_N - dG_(N-1) - dG_1) / RT )
    ! with dG_N "is the change in free enthalpy in the reaction of formation
    ! of 1 mol of clusters of size N from N mol of monomers."
    ! - Gail & Sedlmayr 2013, sec. 13.4.1
    ! See Clouet 2010 https://arxiv.org/abs/1001.4131v2
    ! pages 15+ for a more detailed derivation of
    ! the Gibbs free enegery of the system.
    ! This assumes the clusters to be dilute compared to the total gas
    ! which is resonable

    !! UNCORRECTED
    ! ngas = everything besides the clusters, but as clusters are assumed to be dilute
    ! their number density can be neglected compared to the total
    ! ngas = sum(n(1:nmols))
    ! if(monomer_idx == idx_TiO2 )then
    !   ngas = sum(n(1:nmols))-sum(n(idx_TiO2:idx_Ti10O20))
    ! else
    !   print *, "Clusters other than TiO2 are not yet defined"
    ! endif

    ! print *, "temp", temperature
    gibbs_big = gibbs_free_energy(monomer_idx, cluster_size, temperature) ! kJ*mol**(-1)
    gibbs_small = gibbs_free_energy(monomer_idx, cluster_size-1, temperature)! kJ*mol**(-1)
    gibbs_monomer = gibbs_free_energy(monomer_idx, 1, temperature)! kJ*mol**(-1)
    ! print *, gibbs_big, gibbs_small, gibbs_monomer
    ! correction to the Gibbs free enegery under non-standard pressure of 1 bar.
    ! This only differs in the translational partition function.
    ! ! total gas pressure in units of 1 bar
    ! pressure_scaled = ngas * boltzmann_erg * temperature * 1.e-6
    ! gibbs_corr = temperature * Rgas_kJ * log(pressure_scaled)
    ! gibss correction needs to be added to each gibss energy but _big and _small cancel
    ! The gibbs_corr factor ultimately cancels out ngas and reduced to:
    ! gibbs_corr = temperature * Rgas_kJ *log(pressure_scaled)
    non_standard_correction = (1.e-6_dp * boltzmann_erg * temperature)**(-1)

    gibbs_part = exp( (gibbs_big - gibbs_small - gibbs_monomer)&
                / ( Rgas_kJ * temperature ) )

    k_growth = cluster_growth_rate(monomer_idx, cluster_size-1, temperature)


    !! UNCORRECTED rate
    ! rate = k_growth * ngas * gibbs_part

    ! corrected rate
    rate = k_growth * gibbs_part * non_standard_correction! s^(-1)

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
    integer :: molecule_idx

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
    !    print *, "TIO",cluster_size,gibbs, T, temperature
    else
        if (cluster_size == 1) then
            molecule_idx = monomer_idx
        else
            if(monomer_idx == idx_SiO) then
                molecule_idx = idx_SI2O2 + cluster_size -2
            else if (monomer_idx == idx_MgO) then
                molecule_idx = idx_Mg2O2 + cluster_size -2
            else if (monomer_idx == idx_Al2O3) then
                molecule_idx = idx_Al4O6 + cluster_size -2
            end if
        endif

      gibbs = revHS(T,molecule_idx)*T*Rgas_kJ
    !   print *, "revHS",gibbs,monomer_idx,molecule_idx, T

        ! print *, "There is no thermochemical data on &
        ! clusters ofther than TiO2."
    end if

  end function gibbs_free_energy
#ENDIFKROME
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


  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
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

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  ! k_rev = k_for * revKc
  ! Note that reaction constant revKc is calculated with
  ! reactants and products from reverse reaction
  function revKc(Tgas,ridx,pidx)
    use krome_constants
    use krome_commons
    implicit none
    real*8::revKc,Tgas,dgibss,stoichiometricChange
    integer::ridx(:),pidx(:),i

    ! when considering forward reaction:
    ! Kc = (P°)**(p+p-r-r) * exp(-dGibss_forward°)
    ! where ° means at standard conditions of
    ! P° = 1 bar = (kb*T/1e6) dyn/cm^2 (cgs)
    ! when considering reverse:
    ! 1/Kc = revKc = (kb*T/1e6)**(p+p-r-r) * exp(-dGibss_reverse°)
    ! kb*T/1e6 is to go from 1 atm pressure to number density cm^-3
    ! When not at standard pressure this does not change:
    ! revKc = P**(p+p-r-r) *exp(-dGibss_reverse° - (p+p-r-r)*ln(P/P°))
    !       = (P°)**(p+p-r-r) * exp(-dGibss_reverse°)

    dgibss = 0.d0 ! Gibbs free energy/(R*T)
    stoichiometricChange = 0d0

    do i=1,size(pidx)
       dgibss = dgibss + revHS(Tgas,pidx(i))
       stoichiometricChange = stoichiometricChange + 1
    end do

    do i=1,size(ridx)
       dgibss = dgibss - revHS(Tgas,ridx(i))
       stoichiometricChange = stoichiometricChange - 1
    end do

     revKc = (boltzmann_erg * Tgas * 1e-6)**(stoichiometricChange)&
         * exp(-dgibss)

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    use krome_constants
    use krome_fit
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::Tnist,Tnist2,Tnist3,Tnist4,invTnist,invTnist2,lnTnist
    real*8::xtable(200),multable
#KROME_var_reverse
    integer::idx

    p(:) = 0.d0
    p1_nasa(:,:) = 0.d0
    p2_nasa(:,:) = 0.d0
    Tlim_nasa(:,:) = 0.d0
    p1_nist(:,:) = 0.d0
    p2_nist(:,:) = 0.d0
    Tlim_nist(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)
    ! NIST polynomials are quite differernt
    ! it doesn't like easy stuff...
    Tnist = Tgas * 1.d-3
    Tnist2 = Tnist * Tnist
    Tnist3 = Tnist2 * Tnist
    Tnist4 = Tnist3 * Tnist2
    invTnist = 1d0/Tnist
    invTnist2 = invTnist * invTnist
    lnTnist = log(Tnist)

    hasThermoTable(:) = .false.
#IFKROME_useThermoTables
    yThermoTable(:,:) = 0.d0
    xtable(:) = thermo_tab_Tgas
    multable = thermo_mult_Tgas
#KROME_thermo_tables
#ENDIFKROME

#KROME_kc_reverse_nasa
#KROME_kc_reverse_nist

    !use thermochemical table if available
    if (hasThermoTable(idx)) then
      !NOTE: currently not extrapolation outside tgas table limits
      revHS = fit_anytab1D(xtable, yThermoTable(idx,:), multable, Tgas)
      revHS = revHS/(Rgas_kJ * Tgas)

    ! pick NASA data if present for species
    else if (Tlim_nasa(idx,2) /= 0.d0) then
      !select set of NASA polynomials using temperature
      if(Tlim_nasa(idx,1).le.Tgas .and. Tgas.le.Tlim_nasa(idx,2)) then
         p(:) = p1_nasa(idx,:)

      else if(Tlim_nasa(idx,2)<Tgas .and. Tgas.le.Tlim_nasa(idx,3)) then
         p(:) = p2_nasa(idx,:)

      ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NASA polynomials for enthalpy and enthropy (unitless)
      H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
           p(5)*Tgas4*0.2d0 + p(6)*invT
      S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
           p(5)*Tgas4*0.25d0 + p(7)

      revHS = H - S

    ! else pick NIST data (if present)
    else if (Tlim_nist(idx,2) /= 0.d0) then
      if (Tlim_nist(idx,1) < Tgas .and. Tgas < Tlim_nist(idx,2)) then
        p(:) = p1_nist(idx,:)

      else if (Tlim_nist(idx,2) < Tgas .and. Tgas < Tlim_nist(idx,3)) then
        p(:) = p2_nist(idx,:)

      ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NIST polynomials for enthalpy and enthropy
      ! H in (kJ/mol)
      H = p(1)*Tnist + p(2)*0.5d0*Tnist2 + p(3)*Tnist3/3.d0 + p(4)*Tnist4*0.25d0&
           - p(5)*invTnist + p(6)
      !  Unitsless
      H = H / (Rgas_kJ * Tgas)

      ! S in (J/mol*K)
      S = p(1)*lnTnist + p(2)*Tnist + p(3)*Tnist2*0.5d0 + p(4)*Tnist3/3.d0&
           - p(5)*invTnist2*0.5d0 + p(7)
      !  Unitless. Note: do not use Tnist
      S = S / Rgas_J

      revHS = H - S

    ! return zero is no data exists
    else
      print *, "No thermochemical data of species index", idx
      revHS = 0.d0

    end if


  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
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
    use krome_getphys
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
    use krome_getphys
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

#IFKROME_useH2dust_constant
  !********************************
  !H2 formation on dust using Jura rate
  !dust2gas_ratio in terms of D_solar
  !Usually D/D_sol = Z/Z_sol
  function H2_dustJura(n)
    use krome_commons
    use krome_getphys
    use krome_user_commons
    implicit none
    real*8::n(nspec),H2_dustJura
    real*8::ntot

    ntot = get_Hnuclei(n(:))

    H2_dustJura = n(idx_H)*ntot*3.5d-17*dust2gas_ratio*clump_factor

  end function H2_dustJura
#ENDIFKROME

#IFKROME_useLAPACK
  !*********************************
  subroutine mydgesv(n,Ain,Bin, parent_name)
    !driver for LAPACK dgesv
    integer::n,info,i,ipiv(n)
    real*8,allocatable::tmp(:)
    real*8::A(n,n),B(n),Ain(:,:),Bin(:),suml,sumr,tmpn(n)
    character(len=*)::parent_name
    A(:,:) = Ain(1:n,1:n)
    B(:) = Bin(1:n)
    call dgesv(n,1,A,n,ipiv,B,n,info)
    Bin(1:n) = B(:)

    !write some info about the error and stop
    if(info > 0) then
       allocate(tmp(size(Bin)))
       print *,"ERROR: matrix exactly singular, U(i,i) where i=",info
       print *,' (called by "'//trim(parent_name)//'" function)'

       !dump the input matrix to a file
       open(97,file="ERROR_dump_dgesv.dat",status="replace")
       !dump size of the problem
       write(97,*) "size of the problem:",n
       write(97,*)

       !dump matrix A
       write(97,*) "Input matrix A line by line:"
       do i=1,size(Ain,1)
          tmp(:) = Ain(i,:)
          write(97,'(I5,999E17.8e3)') i,tmp(:)
       end do

       !dump matrix A
       write(97,*)
       write(97,*) "Workin matrix A line by line:"
       do i=1,n
          tmpn(:) = Ain(i,1:n)
          write(97,'(I5,999E17.8e3)') i,tmpn(:)
       end do

       !dump matrix B
       write(97,*)
       write(97,*) "Input/output vector B element by element"
       do i=1,n
          write(97,*) i, Bin(i),B(i)
       end do

       !dump info on matrix A rows
       write(97,*)
       write(97,*) "Info on matrix A rows"
       write(97,'(a5,99a17)') "idx","minval","maxval"
       do i=1,size(Ain,1)
          write(97,'(I5,999E17.8e3)') i, minval(Ain(i,:)), &
               maxval(Ain(i,:))
       end do

       !dump info on matrix sum left and right
       write(97,*)
       write(97,*) "Info on matrix A, sum left/right"
       write(97,'(a5,99a17)') "idx","left","right"
       suml = 0d0
       sumr = 0d0
       do i=1,size(Ain,1)
          if(i>1) suml = sum(Ain(i,:i-1))
          if(i<n) sumr = sum(Ain(i,i+1:))
          write(97,'(I5,999E17.8e3)') i, suml, sumr
       end do
       close(97)

       print *,"Input A and B dumped in ERROR_dump_dgesv.dat"

       stop
    end if

    !if error print some info and stop
    if(info<0) then
       print *,"ERROR: input error position ",info
       print *,' (called by "'//trim(parent_name)//'" function)'
       stop
    end if

  end subroutine mydgesv
#ENDIFKROME

  ! ************************************
  ! solves linear least squares
  subroutine llsq(n, x, y, a, b)

    !****************************************************
    !
    !! LLSQ solves a linear least squares problem matching a line to data.
    !
    !  Discussion:
    !
    !    A formula for a line of the form Y = A * X + B is sought, which
    !    will minimize the root-mean-square error to N data points
    !    ( X(I), Y(I) );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    In: N, the number of data values.
    !
    !    In: X(N), Y(N), the coordinates of the data points.
    !
    !    Out: A, B, the slope and Y-intercept of the
    !    least-squares approximant to the data.
    !
    implicit none
    integer,intent(in)::n
    real*8,intent(out)::a, b
    real*8,intent(in)::x(n), y(n)
    real*8::bot, top, xbar, ybar

    ! special case
    if(n == 1) then
       a = 0d0
       b = y(1)
       return
    end if

    ! average X and Y
    xbar = sum(x) / n
    ybar = sum(y) / n

    ! compute beta
    top = dot_product(x(:) - xbar, y(:) - ybar)
    bot = dot_product(x(:) - xbar, x(:) - xbar)

    ! if top is zero a is zero
    if(top==0d0) then
       a = 0d0
    else
       a = top / bot
    end if

    b = ybar - a * xbar

  end subroutine llsq

end module krome_subs
