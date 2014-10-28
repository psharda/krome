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

    print *,"Initializing dust..."

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
  !load the Qabs from the file fname
  ! for the dust type jtype-th
  subroutine dust_load_Qabs(fname,jtype)
    use krome_commons
    use krome_constants
    implicit none
    character(*)::fname
    integer::ios,nE,i,jtype,nd,nlow,nup
    integer,parameter::nEmax=int(1e4)
    real*8::rout(5),E,asize,Qabs,asize_old,Qabs_old
    real*8::Qabs_tmp(ndust,nEmax),Qabs_E_tmp(nEmax),asize0

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
       Qabs_E_tmp(nE) = E !store energy found
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

    !if not allocated allocate
    if(.not.(allocated(dust_Qabs))) then
       !allocate common Qabs
       allocate(dust_Qabs(ndust,nE))
       allocate(dust_Qabs_E(nE))
       allocate(dust_intBB(ndust,dust_nT))
       allocate(dust_intBB_dT(ndust,dust_nT))
       allocate(dust_intBB_sigma(ndust,dust_nT))
       dust_intBB_Tbb(:) = 0d0 !default
    end if

    !store the number of energies in the common
    dust_Qabs_nE = nE

    !copy temp Qabs to common Qabs
    nlow = (jtype-1)*nd+1
    nup = jtype*nd
    dust_Qabs(nlow:nup,:) = Qabs_tmp(nlow:nup,1:nE)
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
    real*8::Tbb,E,intBB,dE,intBB_dT

    !if already initialized no need to reload
    if(dust_intBB_Tbb(1)/=dust_intBB_Tbb(dust_nT)) return

    !loop on dust bins
    do k=1,ndust
       !loop on Tbb
       do i=1,dust_nT
          !increase Tbb (TbbMax is log and common)
          !Tbb = 1d1**((i-1)*TbbMax/(dust_nT-1))
          Tbb = (i-1) / TbbMult  + TbbMin
          dust_intBB_Tbb(i) = Tbb !store Tbb
          intBB = 0d0
          intBB_dT = 0d0
          !integral Q(E,a)*B(E,Tbb) and Q*dB/dTdust
          do j=2,dust_Qabs_nE
             dE = dust_Qabs_E(j) - dust_Qabs_E(j-1)
             intBB = intBB + 0.5d0 * dE &
                  * (planckBB(dust_Qabs_E(j),Tbb) &
                  * dust_Qabs(k,j) + planckBB(dust_Qabs_E(j-1),Tbb) &
                  * dust_Qabs(k,j-1))
             intBB_dT = intBB_dT + 0.5d0 * dE &
                  * (planckBB_dT(dust_Qabs_E(j),Tbb) &
                  * dust_Qabs(k,j) + planckBB_dT(dust_Qabs_E(j-1),Tbb) &
                  * dust_Qabs(k,j-1))
          end do
          !integral / hplanck (erg/s/cm2)
          dust_intBB(k,i) = intBB * pi * iplanck_eV * eV_to_erg
          !integral / hplanck (erg/s/cm2/K)
          dust_intBB_dT(k,i) = intBB_dT * pi * iplanck_eV * eV_to_erg
          !normalized integral
          dust_intBB_sigma(k,i) = pi * dust_intBB(k,i) &
               / (stefboltz_erg*Tbb**4)
       end do
    end do

  end subroutine dust_init_intBB

  !*************************
  !return the Qabs for the dust bin jdust
  ! and for the given energy in eV
  function get_Qabs(energy,jdust)
    use krome_commons
    implicit none
    integer::jdust,j
    real*8::get_Qabs,energy

    !loop to find the energy and interpolate
    do j=2,dust_Qabs_nE
       if(energy<dust_Qabs_E(j)) then
          get_Qabs = (energy - dust_Qabs_E(j-1)) &
               / (dust_Qabs_E(j) - dust_Qabs_E(j-1)) &
               * (dust_Qabs(jdust,j) - dust_Qabs(jdust,j-1)) &
               + dust_Qabs(jdust,j-1)
          !when found return
          return
       end if
    end do

    !raise an error if nothing found
    print *,"ERROR: no value found for get_Qabs"
    print *,"energy (eV):", energy
    print *,"jdust:",jdust
    stop

  end function get_Qabs

#IFKROME_usePhotoDust
  !***********************
  !compute the integral J(E)*Qabs(E) over the photobins
  ! for the dust bin jdust
  function get_int_JQabs(jdust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_int_JQabs,intJ
    integer::i,j,jdust

    !loop over photo bins
    intJ = 0d0
    do i=1,nPhotoBins
       intJ = intJ + photobinJ(i) * get_Qabs(photoBinEmid(i),jdust) &
            * photoBinEdelta(i)
    end do

    !returns erg/cm2/s
    get_int_JQabs = intJ * iplanck_eV * eV_to_erg * pi

  end function get_int_JQabs
#ENDIFKROME_usePhotoDust

  !*******************************
  !returns the integral in erg/s/cm2 of the BB at
  ! temperature Tbb multiplied by Qabs for the
  ! dust bin jdust
  function get_dust_intBB(jdust,Tbb)
    use krome_commons
    implicit none
    integer::jdust,ibb
    real*8::Tbb,get_dust_intBB

    ibb = (Tbb-TbbMin) * TbbMult + 1

    get_dust_intBB = (Tbb - dust_intBB_Tbb(ibb)) &
         / (dust_intBB_Tbb(ibb+1) - dust_intBB_Tbb(ibb)) &
         * (dust_intBB(jdust,ibb+1) - dust_intBB(jdust,ibb)) &
         +  dust_intBB(jdust,ibb)

  end function get_dust_intBB

  !*******************************
  !returns the integral in erg/s/cm2/K of the dBB/dTdust
  ! at temperature Tbb multiplied by Qabs for the
  ! dust bin jdust
  function get_dust_intBB_dT(jdust,Tbb)
    use krome_commons
    implicit none
    integer::jdust,ibb
    real*8::Tbb,get_dust_intBB_dT

    ibb = (Tbb-TbbMin) * TbbMult + 1

    get_dust_intBB_dT = (Tbb - dust_intBB_Tbb(ibb)) &
         / (dust_intBB_Tbb(ibb+1) - dust_intBB_Tbb(ibb)) &
         * (dust_intBB_dT(jdust,ibb+1) - dust_intBB_dT(jdust,ibb)) &
         +  dust_intBB_dT(jdust,ibb)

  end function get_dust_intBB_dT

  !*******************************
  !compute beta escape using the planck opacity
  function besc(n,Tgas,lj,rhogas)
    use krome_commons
    use krome_subs
    implicit none
    real*8::n(:),besc,Tgas,tau,lj,tau_d,tau_g,rhogas,Tff
    integer::j

    besc = 1d0
    tau_d = 0d0
    do j=1,ndust
       Tff = krome_dust_T(j)
#IFKROME_usedTdust
       Tff = n(nmols+ndust+j)
#ENDIFKROME_usedTdust
       !when temperature difference is small
       ! Tgas can be used instead of Tdust
       if(abs(Tff-Tgas)<5d0) Tff = Tgas
       tau_d = tau_d + xdust(j) * kpla_dust(Tff,j)
    end do

    tau_g = 0d0
    !tau_g = rhogas * 1d1**fit_anytab2D(mayer_x(:), mayer_y(:), mayer_z(:,:), &
    !     mayer_xmul, mayer_ymul,log10(rhogas),log10(Tgas))
    tau = (tau_d + tau_g) * lj

    if(tau<1d0) return
    besc = tau**-2

  end function besc

  !********************************
  !planck opacity for the dust in the
  ! bin jdust for a given temperature Tdust
  function kpla_dust(Tdust,jdust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::kpla_dust,Tdust,intBB
    integer::jdust,ibb

    !ibb = (dust_nT - 1) * log10(Tdust) / TbbMax + 1
    ibb = (Tdust-TbbMin) * TbbMult + 1
    !loop on the temperatures to find the corresponding interval

    intBB = (Tdust-dust_intBB_Tbb(ibb)) &
         / (dust_intBB_Tbb(ibb+1)-dust_intBB_Tbb(ibb)) &
         * (dust_intBB_sigma(jdust,ibb+1)-dust_intBB_sigma(jdust,ibb)) &
         + dust_intBB_sigma(jdust,ibb)

    kpla_dust = intBB * krome_dust_asize2(jdust)

  end function kpla_dust

#IFKROME_usedTdust
  !***********************************
  !dust temperture differential
  function get_dTdust(n,dTgas,vgas,ntot)
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    real*8::get_dTdust(ndust),n(:),dTgas,intBB,pre
    real*8::vgas,ntot,Tgas,rhogas,ljeans,be,m(nspec)
    real*8::fact
    integer::i

    fact = 0.5d0
    m(:) = get_mass()
    Tgas = n(idx_Tgas)
    pre = 0.5d0*fact*vgas*boltzmann_erg*ntot
    rhogas = sum(n(1:nmols)*m(1:nmols))
    ljeans = get_jeans_length_rho(n(:),Tgas,rhogas)
    be = besc(n(:),Tgas,ljeans,rhogas)

    do i=1,ndust
       intBB = get_dust_intBB_dT(i,n(nmols+ndust+i)) * be
       get_dTdust(i) = pre * dTgas / (pre + intBB)
    end do

  end function get_dTdust
#ENDIFKROME_usedTdust

  !*********************************
  !This subroutine computes the dust temperature for each bin
  ! and copies the cooling in the common variable
  subroutine compute_Tdust(n,Tgas)
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    integer::i,j1,j2,jmid
    real*8::Td1,Td2,fact,vgas,ntot,n(:),be,ljeans,rhogas
    real*8::f1,f2,fmid,pre,Tdmid,Tgas,dustCooling,intCMB
    real*8::m(nspec),intJflux

    !compute dust cooling pre-factor (HM79)
    fact = 0.5d0
    ntot = sum(n(1:nmols))
    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas
    pre = 0.5d0*fact*vgas*boltzmann_erg * ntot
    !jeans length cm
    ljeans = get_jeans_length_rho(n(:),Tgas,rhogas)
    !escape
    be = besc(n(:),Tgas,ljeans,rhogas)

    !init dust cooling
    dustCooling = 0d0 

    !init external radiation flux
    intJflux = 0d0

    !loop on dust bins
    do i=1,ndust
       j1 = 1 !first index
       j2 = dust_nT !last index
       !no need to compute Tdust when small amount of dust
       if(n(nmols+1)<1d-18*ntot) then
#IFKROME_usedTdust
          n(nmols+ndust+i) = Tgas
          cycle
#ENDIFKROME_usedTdust
          krome_dust_T(i) = Tgas
          cycle
       end if
       !get initial temperatures
       Td1 = dust_intBB_Tbb(j1)
       Td2 = dust_intBB_Tbb(j2)-1d0
       !compute Tcmb
       intCMB = get_dust_intBB(i,phys_Tcmb)
#IFKROME_usePhotoDust
       !compute external radiation term
       intJflux = get_int_JQabs(i)
#ENDIFKROME_usePhotoDust
       !bisection method
       do 

          !f(x) evaluated at j1 and j2
          f1 = (get_dust_intBB(i,Td1) - intCMB - intJflux) * be &
               - pre * (Tgas - Td1)
          f2 = (get_dust_intBB(i,Td2) - intCMB - intJflux) * be &
               - pre * (Tgas - Td2)

          !compute Tdmid
          Tdmid = .5d0 * (Td1 + Td2)
#IFKROME_usedTdust
          n(nmols+ndust+i) = Tdmid
#ENDIFKROME_usedTdust
          krome_dust_T(i) = Tdmid
          fmid = (get_dust_intBB(i,Tdmid) - intCMB - intJflux) * be &
               - pre * (Tgas - Tdmid)

          !check signs and assign Tdmid
          if(f1*fmid<0d0) then
             Td2 = Tdmid
          else
             Td1 = Tdmid
          end if
          !convergence criterium
          if(abs(Td1-Td2)<1d-8) exit
       end do

#IFKROME_usedTdust
       !no need to compute cooling if dTdust/dT enabled
       cycle
#ENDIFKROME_usedTdust

       !compute the cooling (avoid the difference Tgas-Tdust)
       dustCooling = dustCooling + (get_dust_intBB(i,krome_dust_T(i)) &
            - intCMB - intJflux) * be * n(nmols+i) * krome_dust_asize2(i)
    end do

#IFKROME_usedTdust
    !no need to compute cooling if dTdust/dT enabled
    return
#ENDIFKROME_usedTdust

    !copy (isotropic) cooling
    dust_cooling = 4d0 * pi * dustCooling

  end subroutine compute_Tdust

  !*****************************
  !computes the dust cooling using the temperature difference
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

#IFKROME_useChemisorption
  !load the chemisorption rates from file
  !************************
  subroutine init_chemisorption_rates()
    use krome_commons
    implicit none
    real*8::rout
    integer::ios,icount
    character(len=80)::cout
    character(len=2)::mode
    logical::allocated
    allocated = .false. !determine if tables are already allocated
    mode = "--" !default mode
    print *,"reading chemisorption rates..."
    !open rate file
    open(34,file="surface_chemisorption_rates.dat",status="old")
    do
       read(34,*,iostat=ios) cout !read line as a string
       if(ios==-1) exit !on EOF break loop
       read(cout,*,iostat=ios) rout !convert to double to check if string
       !on error evaluate string
       if(ios.ne.0) then
          cout = trim(cout) !trim line
          !look for non-comment string (i.e. mode)
          if(cout(1:1)/="#") then
             mode = cout(1:2) !store mode
             icount = 0 !starts to count lines
          end if
          cycle !go to the next line
       end if
       icount = icount + 1 !count lines
       !first line is the number of Tdust steps, hence allocate arrays
       if(icount==1.and..not.allocated) then
          dust_rateChem_xsteps = rout
          allocate(dust_rateChem_PC(dust_rateChem_xsteps))
          allocate(dust_rateChem_CP(dust_rateChem_xsteps))
          allocate(dust_rateChem_CC(dust_rateChem_xsteps))
          allocate(dust_rateChem_x(dust_rateChem_xsteps))
          allocated = .true.
       end if
       !second line is the minimum Tdust
       if(icount==2) dust_rateChem_xmin = rout
       !third line is the step in Tdust (linear)
       if(icount==3) then
          dust_rateChem_dx = rout
          dust_rateChem_invdx = 1d0/rout
          dust_rateChem_xfact = (dust_rateChem_xsteps-1) / dust_rateChem_dx &
               / dust_rateChem_xsteps
       end if
       !other lines are data
       if(icount>3) then
          !different modes are differents processes
          if(mode=="PC") then
             dust_rateChem_PC(icount-3) = rout
          elseif(mode=="CP") then
             dust_rateChem_CP(icount-3) = rout
          elseif(mode=="CC") then
             dust_rateChem_CC(icount-3) = rout
          else
             if(icount==4) print *,"WARNING: chemisorption mode "//mode//" skipped!"
          end if
          !store the dependent variable (Tdust)
          dust_rateChem_x(icount-3) = dust_rateChem_xmin + (icount-3-1) * dust_rateChem_dx
       end if
    end do
    close(34)

    print *,"done!"

  end subroutine init_chemisorption_rates
#ENDIFKROME_useChemisorption

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

  !******************
  !dust evaporation following Stahler+1981, eqn.5
  ! ebind: binding energy in K
  function dust_evap(adust,nndust,Tdust,ebind)
    implicit none
    real*8::dust_evap,adust,nndust,Tdust,ebind

    !l0*nu0=1.8d5 cm/s from Stahler+1981
    dust_evap = 3d0 * 1.8d5 * exp(-ebind/Tdust) / adust * nndust

  end function dust_evap

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
