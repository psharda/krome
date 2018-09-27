module krome_dust
contains

#IFKROME_useDust
  subroutine set_dust_distribution(ngas,dust_gas_ratio,alow,aup,phi)
    !krome_init_dust: initialize the dust ditribution
    ! and the dust bin mean sizes. Arguments are
    ! alow the size of the smallest
    ! bin size, aup the largest, phi the exponent
    ! of the MRN power law. These parameters
    ! are optional and can be omitted during the call.
    use krome_commons
    use krome_subs
    use krome_constants
    use krome_getphys
    implicit none
    real*8::rhogas,dmass,ngas(nmols),dust_gas_ratio,d2g
    real*8::iphi1,c,phi1,abin(ndust+1),mass(nspec),myc,phi4
    real*8::alow,aup,phi,dtg(ndustTypes),adust(ndust),Tbb,myx,a0,a1
    real*8::alogmin,alogmax,keyFrac(ndustTypes),frac(ndustTypes)
    integer::i,j,ilow,iup,imax,nd

    d2g = dust_gas_ratio

#KROME_dustPartnerIndex

    phi1 = phi + 1.d0
    phi4 = phi + 4.d0
    iphi1 = 1.d0/phi1

    alogmin = log10(alow)
    alogmax = log10(aup)

    mass(:) = get_mass()
    nd = ndust/ndustTypes

    !dust grain density (g/cm3)
#KROME_dust_grain_density

    !dust binding energy (K)
#KROME_dust_binding_energy

    !key element fraction
#KROME_dust_key_fraction

    frac(:) = keyFrac(:)/sum(keyFrac)

    rhogas = sum(mass(1:nmols)*ngas(:)) !total gas density g/cm3

    !loop on dust types
    do j=1,ndustTypes
       !integration constant
       myc = frac(j)*rhogas*d2g/(4d0*pi/3d0*krome_grain_rho(j) &
            *(aup**phi4-alow**phi4)/phi4)
       ilow = nd * (j - 1) + 1 !lower index
       iup = nd * j !upper index

       !loop to find limits
       do i=1,nd+1
          abin(i) = 1d1**((i-1)*(alogmax-alogmin)/nd+alogmin)
       end do

       !loop to find mean size, span, and dust number density
       do i=1,nd
          adust(i+ilow-1) = ((abin(i+1)**phi4-abin(i)**phi4)/(abin(i+1)**phi1-abin(i)**phi1)*(phi1/phi4))**(1./3.) !(abin(i)+abin(i+1)) * 0.5d0 !mean bin size
          krome_dust_aspan(i+ilow-1) = abin(i+1) - abin(i) !bin span
          xdust(i+ilow-1) = myc*(abin(i+1)**phi1-abin(i)**phi1)/phi1
       end do
#IFKROME_useDustEvol
       !evaluate dust-parnter ratio (e.g. 1dust=1e2 C atoms)
       krome_dust_partner_ratio(ilow:iup) = adust(ilow:iup)**3 &
            * krome_grain_rho(j) / mass(krome_dust_partner_idx(j))
       krome_dust_partner_ratio_inv(ilow:iup) = 1d0 &
            / krome_dust_partner_ratio(ilow:iup)
#ENDIFKROME_useDustEvol
    end do

    krome_dust_asize(:) = adust(:) !store mean size
    krome_dust_asize2(:) = adust(:)**2 !store mean size square
    krome_dust_asize3(:) = adust(:)**3 !store mean size cube

#IFKROME_useDustEvol
    !compute the mass of the dust partner
    do j=1,ndustTypes
       krome_dust_partner_mass(j) =  mass(krome_dust_partner_idx(j))
    end do
#ENDIFKROME_useDustEvol

    !default dust temperature
    krome_dust_T(:) = 3d1

    !init optical properties
#KROME_init_Qabs

    !init integral Qabs(nu)*B(nu)*dnu
#KROME_opt_integral

  end subroutine set_dust_distribution

  !****************************
  !load the Qabs from the file fname
  ! for the dust type jtype-th
  subroutine dust_load_Qabs(fname,jtype)
    use krome_commons
    use krome_constants
    implicit none
    character(*)::fname
    integer::ios,nE,i,jtype,nd,nlow,nup,na,idx,j
    integer,parameter::nEmax=int(1e4),namax=int(1e4)
    real*8::rout(5),E,asize,Qabs,asize_old,Qabs_old
    real*8::Qabs_tmp(namax,nEmax),Qabs_E_tmp(nEmax)
    real*8::Qabs_a_tmp(namax)

    !number of dust bins per type
    nd = ndust/ndustTypes

    !check if Qabs for the given type is already filled with data
    if(Qabs_allocated(jtype)) return
    Qabs_allocated(jtype) = .true.

    !set initial values and counters
    nE = 0
    na = 0
    asize_old = -1d99
    Qabs_old = -1d99

    !open file and read
    open(32,file=fname,status="old",iostat=ios)
    !loop on file lines
    do
       read(32,*,iostat=ios) rout(:)
       if(ios<0) exit !eof
       !skip when # is found
       if(ios==59.or.ios==5010) cycle
       !skip blanks
       if(ios/=0) cycle

       !read values
       E = rout(2) * planck_eV !eV
       Qabs = rout(3) !dimensionless
       asize = rout(1) !cm

       !check for next energy block
       if(asize/=asize_old) then
          nE = 0 !restart energy counter
          na = na + 1 !increase size counter
       end if

       nE = nE + 1 !increase energy index
       Qabs_E_tmp(nE) = E !store energy found, eV
       Qabs_a_tmp(na) = asize !store size found, cm
       Qabs_tmp(na,nE) = Qabs !store Qabs
       asize_old = asize
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

    !upper and lower limits
    nlow = (jtype-1)*nd+1
    nup = jtype*nd
    !copy temp Qabs to common Qabs and interpolate over grain sizes
    !loop on dust bins
    do i=nlow,nup
       idx = -1 !deafult error-rising index
       !loop on opt table sizes to found dust size corresponding index
       do j=1,na
          if(krome_dust_asize(i)<Qabs_a_tmp(j)) then
             idx = j !found index
             exit
          end if
       end do
       !interpolate on the dust size
       dust_Qabs(i,:) = (Qabs_tmp(idx,1:nE) - Qabs_tmp(idx-1,1:nE)) &
            *(krome_dust_asize(i)-Qabs_a_tmp(idx-1)) &
            / (Qabs_a_tmp(idx)-Qabs_a_tmp(idx-1)) + Qabs_tmp(idx-1,1:nE)
    end do

    !copy energies
    dust_Qabs_E(:) = Qabs_E_tmp(1:nE)

  end subroutine dust_load_Qabs

  !***********************
  !This sub prepares a table of the BB integral
  ! weighted by the Qabs as a function of Tbb e asize
  subroutine dust_init_intBB()
    use krome_subs
    use krome_commons
    use krome_constants
    use krome_phfuncs
    implicit none
    integer::i,j,k,ios,nread
    real*8::Tbb,E,intBB,dE,intBB_dT,rout(4)
    logical::exists

    !if already initialized no need to reload
    if(dust_intBB_Tbb(1)/=dust_intBB_Tbb(dust_nT)) return


    !check if tables are already computed
    inquire(file="KROME_dust_intBB.dat", exist=exists)

    !if tables are already present load from file
    if(exists) then
       print *,"Loading dust tables from file..."
       open(33,file="KROME_dust_intBB.dat",status="old")
       nread = 0
       do
          read(33,*,iostat=ios) k,i,rout(:)
          if(ios.ne.0) exit
          nread = nread + 1
          dust_intBB_Tbb(i) = rout(1)
          dust_intBB(k,i) = rout(2)
          dust_intBB_dT(k,i) = rout(3)
          dust_intBB_sigma(k,i) = rout(4)
       end do
       close(33)

       !check if the number of data loaded is OK
       if(nread/=dust_nT*ndust) then
          print *,"ERROR: the size of the file KROME_dust_intBB.dat"
          print *," seems to be wrong. Delete it and restart the "
          print *," executable."
          stop
       end if
       return
    end if

    print *,"Computing dust tables..."
    print *,"dust bins:",ndust
    print *,"(it may take a while)"

    open(33,file="KROME_dust_intBB.dat",status="replace")
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
          dust_intBB_sigma(k,i) = dust_intBB(k,i) &
               / (stefboltz_erg*Tbb**4)
          !write the data on a file
          write(33,'(2I8,99E17.8e3)') k, i, dust_intBB_Tbb(i), &
               dust_intBB(k,i), dust_intBB_dT(k,i), dust_intBB_sigma(k,i)
       end do
    end do
    close(33)

  end subroutine dust_init_intBB

  !************************
  !interpolate dust Qabs over photobins
  subroutine interp_qabs()
    use krome_commons
    implicit none
    integer::i,jdust

#IFKROME_usePhotoDust
    !loop on dust types
    do jdust=1,ndust
       !loop on photobins
       do i=1,nPhotoBins
          !interpolate Qabs on photo bins
          dust_Qabs_interp(i,jdust) = get_Qabs(photoBinEmid(i),jdust)
       end do
    end do
#ENDIFKROME_usePhotoDust

  end subroutine interp_qabs

  !*************************
  !return the Qabs for the dust bin jdust
  ! and for the given energy in eV
  function get_Qabs(energy,jdust)
    use krome_commons
    implicit none
    integer::jdust,j
    real*8::get_Qabs,energy

    !rise error if array size problem
    if(dust_Qabs_nE<=0) then
       print *,"ERROR: when interpolating Qabs found zero energy points!"
       print *,"You should probably intialize dust before initializing"
       print *," photochemistry..."
       stop
    end if

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

    !rise an error if nothing found
    print *,"ERROR: no value found for get_Qabs"
    print *,"energy (eV):", energy
    print *,"energy limits (eV):", dust_Qabs_E(1), dust_Qabs_E(dust_Qabs_nE)
    print *,"jdust:",jdust
    stop

  end function get_Qabs

  !************************
  !jbin is the photobin, jdust the dust bin
  function get_Qabs_bin(jbin,jdust)
    use krome_commons
    implicit none
    real*8::get_Qabs_bin
    integer::jbin,jdust

    get_Qabs_bin = dust_Qabs_interp(jbin,jdust)

  end function get_Qabs_bin

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
       intJ = intJ + photobinJ(i) &
            !* get_Qabs(photoBinEmid(i),jdust) &
            * get_Qabs_bin(i,jdust) &
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
    use krome_phfuncs
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
       if(xdust(j)<1d-20) cycle
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
    besc = tau**(-2)

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

    kpla_dust = intBB * krome_dust_asize2(jdust) * pi

  end function kpla_dust

#IFKROME_usedTdust
  !***********************************
  !dust temperture differential
  function get_dTdust(n,dTgas,vgas,ntot)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_getphys
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
    use krome_getphys
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
       if(xdust(i)<1d-30*ntot) then
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
            - intCMB - intJflux) * be * xdust(i) * krome_dust_asize2(i)
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
  function krome_dust_growth(natom,Tgas,Tdust,vgas,atom_mass,rho0)
    !krome_dust_growth: compute dust growth in cm/s
    implicit none
    real*8::krome_dust_growth,natom,Tgas,Tdust,vgas,atom_mass,rho0

    krome_dust_growth = krome_dust_stick(Tgas,Tdust) * vgas * &
         atom_mass / 4d0 / rho0 * natom

  end function krome_dust_growth

  !*******************
  function krome_dust_stick(Tgas,Tdust)
    !krome_dust_stick: sticking coefficient (Leitch-Devlin & Williams 1985)
    ! (Grassi2012, eqn.26)
    real*8::krome_dust_stick,Tgas,Tdust

    krome_dust_stick = 1.9d-2 * Tgas * (1.7d-3 * Tdust + .4d0) &
         * exp(-7.d-3 * Tgas)
    krome_dust_stick = min(krome_dust_stick,1d0)

  end function krome_dust_stick

  !******************
  !dust evaporation following Evans+93: cm/s
  ! ebind: binding energy in K
  function dust_evap(Tdust,ebind,amass,asize2,rho0)
    use krome_constants
    implicit none
    real*8::dust_evap,Tdust,ebind,nu0,amass
    real*8::asize2,rho0

    dust_evap = 0d0 !default
    nu0 = 1d12 !1/s
    if(asize2<=0d0) return
    dust_evap = nu0 * amass/4d0/pi/asize2/rho0 * exp(-ebind/Tdust)

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
    mgrain = adust**3 *2.3d0 / (p_mass)
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
    !give back the H2 formation efficiency on silicates
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
    !give back the H2 formation efficiency on C
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
    !sticking coefficient for the H2 formation
    real*8::stick,myTgas,myTdust
    stick = 1.d0 / (1.d0 + 0.4*sqrt((myTgas+myTdust)/1d2) &
         + 0.2*myTgas/1d2 + 0.08*(myTgas/1d2)**2)
    if(stick>1.d0 .or. stick<0.d0) then
       print *,"problem on stick coefficient"
       stop
    end if
  end function stick

#ENDIFKROME

#IFKROME_usePhotoDust_3D
  subroutine load_int_JQabs_tab(integrand)
    use krome_commons
    use krome_constants
    !
    real*8,dimension(nPhotoBins) :: integrand ! tGamma in each photobin
    ! table related local variables
    character(len=255) :: fname
    integer :: nrec, lunit, i, j
    real(kind=8) :: energyL, energyR
    real(kind=8), allocatable, dimension(:) :: emid, el, eu, tGamma
    !
    ! open file, get nr of records, and read data
    fname = 'dust_table_absorption.dat'
    nrec = 0
    open(newunit=lunit, file=fname) 
    do 
       read(lunit,*, end=100)
       nrec = nrec + 1
    end do
    100 rewind(lunit)
    allocate(emid(nrec),el(nrec), eu(nrec), tGamma(nrec))
    do j=1,nrec
       read(lunit,*) emid(j), el(j), eu(j), tGamma(j)
    end do
    close(lunit)
    !
    ! use table to find tGamma integrated across each photo bin
    do i=1,nPhotoBins
       energyL = photoBinEleft(i) * eV_to_erg
       energyR = photoBinEright(i) * eV_to_erg !energy of the bin in erg
       integrand(i) = 0.0_8
       do j = 1, nrec
          if (el(j) <= energyR .and. eu(j) >= energyL) then
             if (el(j) >= energyL .and. eu(j) <= energyR) &
               integrand(i) = integrand(i) + tGamma(j)
             if (el(j) >= energyL .and. eu(j) > energyR) &
               integrand(i) = integrand(i) + tGamma(j) * (energyR - el(j)) / (eu(j) - el(j))
             if (el(j) < energyL .and. eu(j) <= energyR) &
               integrand(i) = integrand(i) + tGamma(j) * (eu(j) - energyL) / (eu(j) - el(j))
             if (el(j) < energyL .and. eu(j) > energyR) &
               integrand(i) = integrand(i) + tGamma(j) * (energyR - energyL) / (eu(j) - el(j))
          end if
       end do
    end do
    !
    deallocate(emid,el,eu,tGamma)
  end subroutine load_int_JQabs_tab
  !***********************
  !compute the absorbed radiation by the dust by integrating over the photobins
  !asuming a tabularised dust distribution
  function get_int_JQabs_tab()
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_int_JQabs_tab
    real*8,dimension(nPhotoBins), save :: integrand ! tGamma in each photobin
    logical, save :: first_call = .true.
    !
    ! Load tabularised data on first call.
    ! Make a double-nested threadprivate check of if this is first call
    ! to avoid race-condition without parallel overhead once table has been loaded
    !
    if (first_call) then
       !$omp critical
       if (first_call) call load_int_JQabs_tab(integrand)
       first_call = .true.
       !$omp end critical
    end if

    !loop over photo bins
    get_int_JQabs_tab = sum(photoBinJ * integrand) * eV_to_erg ! make sure result is in erg s^-1 cm^-3 as for table file
    !
  end function get_int_JQabs_tab
  !
  subroutine setup_2D_dust_tables
    use krome_commons
    implicit none
    integer       :: iAv
    real(kind=8)  :: lambda_abs, wl, wu, log10_Av
    logical,      save :: first_call=.true.
    integer,      save :: nAv, nrec
    real(kind=8), save :: log10_Av_lb, dlog10_Av
    real(kind=8), allocatable, dimension(:), save :: Av_tab, lambda_tab
    real(kind=8), allocatable, dimension(:,:,:), save :: dust_tab_H2_3D, dust_tab_cool_3D, dust_tab_Tdust_3D

    ! Compute absorption from radiation field
    lambda_abs = get_int_JQabs_tab()

    ! Make sure 3D tables are loaded
    !
    ! This will also define
    !    dust_tab_ngas, dust_tab_Tgas,
    !    dust_mult_ngas, dust_mult_Tgas
    ! for 2D table operation, and log10_Av_lb, dlog10_Av
    if (first_call) then
       !$omp critical
       if (first_call) call load_tables()
       first_call = .false.
       !$omp end critical
    end if

    ! Translate from absorption to an Av
    log10_Av = log10(convert_to_Av(lambda_abs))

    ! Find index (iAv) and weights(wl,wu; wl+wu=1)
    iAv = floor((log10_Av - log10_Av_lb)/dlog10_Av)+1
    iAv = min(max(iAv,1),nAv)
    wu  = min(max(log10_Av - ((iAv-1)*dlog10_Av + log10_Av_lb),0.0_8),1.0_8)
    wl = 1.0_8 - wu

    ! Interpolate in 3D table to generate 2D tables. Check if we are on edge point for ub (then wl==1.)
    if (iAv < nAv) then
      dust_tab_H2    = dust_tab_H2_3D(:,:,iAv)*wl + dust_tab_H2_3D(:,:,iAv+1)*wu
      dust_tab_cool  = dust_tab_cool_3D(:,:,iAv)*wl + dust_tab_cool_3D(:,:,iAv+1)*wu
      dust_tab_Tdust = dust_tab_Tdust_3D(:,:,iAv)*wl + dust_tab_Tdust_3D(:,:,iAv+1)*wu
    else
      dust_tab_H2    = dust_tab_H2_3D(:,:,nAv)
      dust_tab_cool  = dust_tab_cool_3D(:,:,nAv)
      dust_tab_Tdust = dust_tab_Tdust_3D(:,:,nAv)
    endif
  contains
     subroutine load_tables
       use krome_fit
       implicit none
       integer :: lunit, i, j, nTgas, nngas
       real(kind=8), dimension(:), allocatable :: dust_tab_Av
       character(len=255) :: fname
       ! Load Av-Lambda table
       fname = 'dust_table_Av_Lambda.dat'
       open(newunit=lunit, file=fname)
       nrec = 0
       do
          read(lunit,*, end=200)
          nrec = nrec + 1
       end do
       200 rewind(lunit)
       allocate(Av_tab(nrec), lambda_tab(nrec))
       do j=1,nrec
          read(lunit,*) Av_tab(j), lambda_tab(j)
       end do
       close(lunit)
       ! Load the three 3D tables
       nTgas = 50 ! FIXME, hardcoded!
       nngas = 50
       nAv = 20
       allocate(dust_tab_Tdust_3D(nTgas,nngas,nAv), &
                dust_tab_cool_3D(nTgas,nngas,nAv), &
                dust_tab_H2_3D(nTgas,nngas,nAv), &
                dust_tab_Av(nAv))
       !
       call init_anytab3D("dust_table_cool_3D.dat",dust_tab_ngas(:), &
            dust_tab_Tgas(:), dust_tab_Av(:), &
            dust_tab_cool_3D(:,:,:), dust_mult_ngas, &
            dust_mult_Tgas, dlog10_Av)
       call init_anytab3D("dust_table_Tdust_3D.dat",dust_tab_ngas(:), &
            dust_tab_Tgas(:), dust_tab_Av(:), &
            dust_tab_Tdust_3D(:,:,:), dust_mult_ngas, &
            dust_mult_Tgas, dlog10_Av)
       call init_anytab3D("dust_table_H2_3D.dat",dust_tab_ngas(:), &
            dust_tab_Tgas(:), dust_tab_Av(:), &
            dust_tab_H2_3D(:,:,:), dust_mult_ngas, &
            dust_mult_Tgas, dlog10_Av)

       log10_Av_lb = dust_tab_Av(1)
       dlog10_Av   = dust_tab_Av(2) - dust_tab_Av(1)

       deallocate(dust_tab_Av)
     end subroutine load_tables
     !
     function convert_to_Av(lambda_abs) result(Av)
       implicit none
       real(kind=8), intent(in) :: lambda_abs
       real(kind=8)             :: Av
       !
       real(kind=8) :: w
       integer :: i
       !
       if (lambda_abs >= lambda_tab(1)) then
         Av = Av_tab(1)
         return
       endif
       !
       if (lambda_abs <= lambda_tab(nrec)) then
         Av = Av_tab(nrec)
         return
       endif
       ! linear interpolation for Av in log of lambda_abs
       do i=2,nrec
          if (lambda_abs > lambda_tab(i)) then
            w = (log(lambda_abs) - log(lambda_tab(i))) / (log(lambda_tab(i-1)) - log(lambda_tab(i)))
            Av = Av_tab(i) * (1.0_8 - w) + Av_tab(i-1) * w
            return
          end if
       end do
       allocate(Av_tab(nrec), lambda_tab(nrec))
     end function convert_to_Av
  end subroutine setup_2D_dust_tables
#ENDIFKROME

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_fit
    implicit none

#IFKROME_dust_table_2D
    call init_anytab2D("dust_table_cool.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_cool(:,:), dust_mult_ngas, &
         dust_mult_Tgas)
    call init_anytab2D("dust_table_Tdust.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
         dust_mult_Tgas)
    call init_anytab2D("dust_table_H2.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_H2(:,:), dust_mult_ngas, &
         dust_mult_Tgas)
#ENDIFKROME

#IFKROME_dust_table_3D
    call init_anytab3D("dust_table_cool.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_AvVariable(:), &
         dust_tab_cool(:,:,:), dust_mult_ngas, &
         dust_mult_Tgas, dust_mult_AvVariable)
    call init_anytab3D("dust_table_Tdust.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_AvVariable(:), &
         dust_tab_Tdust(:,:,:), dust_mult_ngas, &
         dust_mult_Tgas, dust_mult_AvVariable)
    call init_anytab3D("dust_table_H2.dat",dust_tab_ngas(:), &
         dust_tab_Tgas(:), dust_tab_AvVariable(:), &
         dust_tab_H2(:,:,:), dust_mult_ngas, &
         dust_mult_Tgas, dust_mult_AvVariable)
#ENDIFKROME

  end subroutine init_dust_tabs

end module krome_dust
