module krome_dust

#IFKROME_useDust
contains

  subroutine set_dust_distribution(alow,aup,phi)
    !krome_init_dust: initialize the dust ditribution
    ! and the dust bin mean sizes. Arguments are
    ! alow the size of the smallest
    ! bin size, aup the largest, phi the exponent
    ! of the MRN power law. These parameters
    ! are optional and can be omitted during the call.
    use krome_commons
    use krome_subs
    implicit none
    real*8::rhogas,dmass,ngas(nmols),n(nspec)
    real*8::iphi1,c,phi1,abin(ndust+1),mass(nspec),myc
    real*8::alow,aup,phi,dtg(ndustTypes),adust(ndust),Tbb,myx,a0,a1
    integer::i,j,ilow,iup,imax,nd

#KROME_dustPartnerIndex

    phi1 = phi + 1.d0
    iphi1 = 1.d0/phi1

    mass(:) = get_mass()
    nd = ndust/ndustTypes
    krome_grain_rho = 2.3d0 !dust grain density g/cm3 (graphite fits all)

    rhogas = sum(mass(1:nmols)*ngas(:)) !total gas density g/cm3

    !WARNING: this algorithm evaluates the bin sizes in order to
    ! have the same amount of dust per each bin, following MNR
    ! distribution. change according to your needs
    do j=1,ndustTypes
       ilow = nd * (j - 1) + 1 !lower index
       iup = nd * j !upper index

       abin(1) = alow !lower limit size
       abin(nd+1) = aup !upper limit
       myc = (aup**phi1-alow**phi1) / phi1 !distribution total integral

       !loop to find limits (analitically)
       do i=2,nd
          abin(i) = (myc*phi1/nd + abin(i-1)**phi1)**iphi1
       end do

       !loop to find mean size 
       do i=1,nd
          adust(i+ilow-1) = (abin(i)+abin(i+1)) * 0.5d0 !mean bin size
          krome_dust_aspan(i+ilow-1) = abin(i+1) - abin(i) !bin span
       end do

       !amount of dust per bin
       xdust(ilow:iup) = 1d0 / nd

       !evaluate dust-parnter ratio (e.g. 1dust=1e2 C atoms)
       krome_dust_partner_ratio(ilow:iup) = adust(ilow:iup)**3 &
            * krome_grain_rho / mass(krome_dust_partner_idx(j))
       krome_dust_partner_ratio_inv(ilow:iup) = 1d0 &
            / krome_dust_partner_ratio(ilow:iup)
    end do

    krome_dust_asize(:) = adust(:) !store mean size
    krome_dust_asize2(:) = adust(:)**2 !store mean size square
    krome_dust_asize3(:) = adust(:)**3 !store mean size cube


    !compute the mass of the dust partner
    do j=1,ndustTypes
       krome_dust_partner_mass(j) =  mass(krome_dust_partner_idx(j))
    end do

    !default dust temperature
    krome_dust_T(:) = 3d1

    !init optical properties
#KROME_init_Qabs

    !init integral Qabs(nu)*B(nu)*dnu
#KROME_opt_integral

    print *,"Dust initialized!"

  end subroutine set_dust_distribution

  !****************************
  subroutine dust_load_Qabs(fname,jtype)
    use krome_commons
    use krome_constants
    implicit none
    character(*)::fname
    integer::ios,nE,i,jtype,nd
    integer,parameter::nEmax=int(1e4)
    real*8::rout(5),E,asize,Qabs,asize_old,Qabs_old
    real*8::Qabs_tmp(ndust,nEmax),Qabs_E_tmp(nEmax),asize0

    !if already loaded quit subroutine
    if(allocated(dust_Qabs)) return
    
    !number of dust bins per type
    nd = ndust/ndustTypes

    !set initial values
    nE = 0
    asize_old = -1d99
    Qabs_old = -1d99
    asize0 = -1d99

    !open file
    open(32,file=fname,status="old",iostat=ios)
    !loop on file lines
    do
       read(32,*,iostat=ios) rout(:)
       if(ios<0) exit !eof
       if(ios/=0) cycle !skip blanks

       !read values
       E = rout(2) * planck_eV !eV
       Qabs = rout(3) !dimensionless
       asize = rout(1) !cm
       nE = nE + 1 !increase energy index
       Qabs_E_tmp(nE) = E
       !loop on dust size to get Qabs
       do i=(jtype-1)*nd+1,jtype*nd
          !interpolate Qabs from asize
          if(krome_dust_asize(i).ge.asize0 .and. &
               krome_dust_asize(i).le.asize) then
             Qabs_tmp(i,nE) = (krome_dust_asize(i)-asize0) &
                  / (asize-asize0) * (Qabs-Qabs_old) + Qabs_old
          end if
       end do

       !check for next energy block
       if(asize/=asize_old) then
          asize0 = asize_old
          nE = 0
       end if
       !store old
       asize_old = asize
       Qabs_old = Qabs
    end do

    close(32)

    !allocate common Qabs
    allocate(dust_Qabs(ndust,nE))
    allocate(dust_Qabs_E(nE))
    allocate(dust_intBB(ndust,dust_nT))

    dust_Qabs_nE = nE

    !copy temp Qabs to common Qabs
    dust_Qabs(:,:) = Qabs_tmp(:,1:nE)
    dust_Qabs_E(:) = Qabs_E_tmp(1:nE)
    
  end subroutine dust_load_Qabs

  !***********************
  !This sub prepares a table of the BB integral
  ! weighted by the Qabs as a function of Tbb e asize
  subroutine dust_init_intBB()
    use krome_subs
    use krome_commons
    use krome_constants
    implicit none
    integer::i,j,k
    real*8::Tbb,E,intBB,dE

    !loop on dust bins
    do k=1,ndust
       !loop on Tbb
       do i=1,dust_nT
          Tbb = 1d1**((i-1)*4d0/(dust_nT-1))
          dust_intBB_Tbb(i) = Tbb
          intBB = 0d0
          !integral Q(E,a)*B(E,Tbb)
          do j=2,dust_Qabs_nE
             dE = dust_Qabs_E(j) - dust_Qabs_E(j-1)
             intBB = intBB + 0.5d0 * dE * (planckBB(dust_Qabs_E(j),Tbb) &
                  * dust_Qabs(k,j) + planckBB(dust_Qabs_E(j-1),Tbb) &
                  * dust_Qabs(k,j-1)) 
          end do
          !integral / hplanck (erg/s/cm2)
          dust_intBB(k,i) = intBB * pi * iplanck_eV * eV_to_erg
       end do
    end do

  end subroutine dust_init_intBB

  !*******************************
  function besc(n,Tgas)
    use krome_commons
    use krome_subs
    implicit none
    real*8::n(:),besc,Tgas,tau
    integer::j

    tau = 0d0
    do j=1,ndust
       tau = tau + xdust(j) * kpla_dust(krome_dust_T(j),j)
    end do
    tau = tau * get_jeans_length(n(:),Tgas)
    besc = min(1d0,(tau+1d-40)**-2)

  end function besc

  !********************************
  function kpla_dust(Tdust,jdust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::kpla_dust,Tdust,intBB
    integer::i,jdust

    do i=1,dust_nT
       if(dust_intBB_Tbb(i)>Tdust) then
          intBB = (Tdust-dust_intBB_Tbb(i-1)) &
               / (dust_intBB_Tbb(i)-dust_intBB_Tbb(i-1)) &
               * (dust_intBB(jdust,i)-dust_intBB(jdust,i-1)) &
               + dust_intBB(jdust,i-1)
          exit
       end if
    end do

    kpla_dust = intBB / (stefboltz_erg*Tdust**4) * pi &
         * krome_dust_asize2(jdust)

  end function kpla_dust
  
  !*********************************
  !This subroutine compute the dust temperature for each bin
  subroutine compute_Tdust(n,Tgas)
    use krome_commons
    use krome_constants
    implicit none
    integer::i,j1,j2,jmid
    real*8::Td1,Td2,fact,vgas,ntot,n(:),be
    real*8::f1,f2,fmid,pre,Tdmid,Tgas,dustCooling

    !compute dust cooling pre-factor (HM79)
    fact = 0.5d0
    ntot = sum(n(1:nmols))
    vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas
    pre = fact*vgas*boltzmann_erg * ntot
    be = besc(n(:),Tgas)

    !init dust cooling
    dustCooling = 0d0 

    !loop on dust bins
    do i=1,ndust
       j1 = 1 !first index
       j2 = dust_nT !last index
       do 
          !compute Tdust in j1 and j2
          Td1 = dust_intBB_Tbb(j1)
          Td2 = dust_intBB_Tbb(j2)
          !f(x) evaluated at j1 and j2
          f1 = dust_intBB(i,j1) * be - pre * (Tgas-Td1)
          f2 = dust_intBB(i,j2) * be - pre * (Tgas-Td2)
          !next j value
          jmid = (j1+j2) / 2
          Tdmid = dust_intBB_Tbb(jmid)
          fmid = dust_intBB(i,jmid) * be - pre * (Tgas-Tdmid)
          !check signs and assign jmid
          if(f1*fmid<0d0) then
             f2 = fmid
             j2 = jmid
          else
             f1 = fmid
             j1 = jmid
          end if
          !when contiguous break loop
          if(abs(j1-j2)<2) exit
       end do
       dustCooling = dustCooling + 4d0 * pi * dust_intBB(i,jmid) * be &
            * xdust(i) * krome_dust_asize2(i)
       !store temperature
       krome_dust_T(i) = .5d0 * (dust_intBB_Tbb(j1) + dust_intBB_Tbb(j2))
    end do

    dust_cooling = dustCooling
    
  end subroutine compute_Tdust
  
  !*****************************
  function dustCool(adust2,nndust,Tgas,Tdust,ntot)
    use krome_constants
    real*8::dustCool,adust2,Tgas,nndust,ntot,fact,vgas,Tdust

    !factor of contribution for species other than protons
    ! mean value, see Hollenbach and McKee 1979 for a
    ! more accurate value
    fact = 0.5d0
    vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas

    dustCool = 2.d0 * boltzmann_erg * nndust * adust2 * &
         fact * vgas * (Tgas - Tdust) * ntot

  end function dustCool

  !******************
  function sgn(arg)
    !return sign of a double arg
    real*8::sgn,arg
    sgn = 1.d0
    if(arg>=0.d0) return
    sgn = -1.d0
  end function sgn


  !*********************
  function krome_dust_grow(nndust,natom,Tgas,Tdust,vgas,adust)
    !krome_dust_grow: compute dust formation in cm3/s 
    ! (Grassi2012, eqn.25)
    use krome_constants
    implicit none
    real*8::krome_dust_grow,nndust,natom,Tgas,Tdust,vgas,adust,seed

    seed = #KROME_dust_seed !1/cm3
    krome_dust_grow = 1d-4 * pi * adust**2 * (max(nndust,0.d0) + seed) &
         * max(natom,0.d0) * krome_dust_stick(Tgas,Tdust) * vgas

  end function krome_dust_grow

  !*******************
  function krome_dust_stick(Tgas,Tdust)
    !krome_dust_stick: sticking coefficient (Leitch-Devlin & Williams 1985)
    ! (Grassi2012, eqn.26)
    real*8::krome_dust_stick,Tgas,Tdust

    krome_dust_stick = 1.9d-2 * Tgas * (1.7d-3 * Tdust + .4d0) &
         * exp(-7.d-3 * Tgas)

  end function krome_dust_stick

  !***************
  function krome_dust_sput(Tgas,adust,natom,nndust)
    use krome_constants
    use krome_commons
    !sputtering rate using nozawa 2006 yields as impact efficiency
    real*8::krome_dust_sput,Tgas,adust,natom,logT,nndust,y,logy
    real*8::a0,a1,a2,a3,mgrain

    if(Tgas<1d5) then
       krome_dust_sput = 0.d0
       return
    end if
    logT = log10(Tgas)
    a0 = -3.9280689440337335d-01
    a1 = 1.9746828032397993d+01
    a2 = -4.0031865960055839d+00
    a3 = 7.8081665612858187d+00

    !yield contains thermal speed (y*vgas)
    logy = exp(-a0 * logT) / (a1 + a2 * logT) + a3
    y = 1d1**logy
    mgrain = adust**3 *krome_grain_rho / (p_mass)
    krome_dust_sput = y*natom*nndust*adust**2 / mgrain

    if(krome_dust_sput>1.d0) then
       print *,"sputtering>1!"
       print *,"sputtering ","adust ","natom ","Tgas ","nndust"
       print *,krome_dust_sput,adust,natom,Tgas,nndust
       stop
    end if

  end function krome_dust_sput

  !*******************************
  function krome_H2_dust(nndust,Tgas,Tdust,nH,H2_eps_f,myvgas)
    !H2 formed on dust (1/cm3/s)
    use krome_constants
    use krome_commons
    real*8::H2_dust, krome_H2_dust,Tgas,Tdust(:)
    real*8::myvgas,H2_eps,nndust(:),nH,H2_eps_f
    integer::i

    H2_dust = 0.d0 
    do i = 1,size(Tdust)
       H2_eps = H2_eps_f(Tgas, Tdust(i))
       H2_dust = H2_dust + 0.5d0 * nH * myvgas * nndust(i) &
            * krome_dust_asize2(i) &
            * pi * H2_eps * stick(Tgas, Tdust(i))
    end do

    krome_H2_dust = H2_dust

  end function krome_H2_dust

  !*************************
  function H2_eps_Si(myTgas, myTdust)
    !restituisce l'epsilon del Si
    real*8::H2_eps_Si,Ec,Es,Ep,apc,func
    real*8::myTgas,myTdust
    Ec = 1.5d4 !K
    Es = -1d3 !K
    Ep = 7d2 !K
    !m (must be in m even if in CS2009 it's in \AA!!, Cazaux2012 private comm.)
    apc = 1.7d-10 
    func = 2.d0 * exp(-(Ep-Es)/(Ep+myTgas)) / (1.d0+sqrt((Ec-Es)/(Ep-Es)))**2
    H2_eps_Si = 1.d0/(1.+16.*myTdust/(Ec-Es) * exp(-Ep/myTdust)&
         * exp(4d9*apc*sqrt(Ep-Es))) + func

    H2_eps_Si = min(H2_eps_Si,1d0)
    if(H2_eps_Si<0.d0) then
       print *,"problem on H2_eps_Si"
       print *,H2_eps_Si
       stop
    end if
  end function H2_eps_Si


  !*************************
  function H2_eps_C(myTgas, myTdust)
    !resituisce l'epsilon per il C
    real*8::H2_eps_C,myTgas,myTdust
    real*8::Ep,Ec,Es,Th
    Ep = 8d2 !K
    Ec = 7d3 !K
    Es = 2d2 !K
    TH = 4.d0*(1.d0+sqrt((Ec-Es)/(Ep-Es)))**(-2) * exp(-(Ep-Es)/(Ep+myTgas))
    H2_eps_C = (1.d0-TH) / (1.d0+0.25*(1.d0+sqrt((Ec-Es)/(Ep-Es)))**2 &
         * exp(-Es/myTdust))
    if(H2_eps_C>1.d0 .or. H2_eps_C<0.d0) then
       print *,"problem on H2_eps_C"
       stop
    end if
  end function H2_eps_C

  !***************************
  function stick(myTgas,myTdust)
    !sticking coefficient per la formazione di H2
    real*8::stick,myTgas,myTdust
    stick = 1.d0 / (1.d0 + 0.4*sqrt((myTgas+myTdust)/1d2) &
         + 0.2*myTgas/1d2 + 0.08*(myTgas/1d2)**2)
    if(stick>1.d0 .or. stick<0.d0) then
       print *,"problem on stick coefficient"
       stop
    end if
  end function stick

#ENDIFKROME
end module krome_dust
