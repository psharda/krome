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


  !***********************
  subroutine dustOptIntegral(dust_opt_Em, dust_opt_Tbb, dust_opt_asize, &
       dust_opt_nu, dust_opt_Qabs)
    use krome_commons
    integer::i,imax,j,nd
    real*8,allocatable::dust_opt_Em(:,:), dust_opt_Tbb(:)
    real*8::dust_opt_asize(:), dust_opt_nu(:), dust_opt_Qabs(:,:),Tbb
    if(allocated(dust_opt_Em)) return !skip evaluation if already done
    nd = size(krome_dust_asize)
    !computes and stores Qabs(nu)*B(nu)*dnu integral for each dust bin
    imax = int(1e4) !number of Tbb values
    allocate(dust_opt_Em(nd,imax), dust_opt_Tbb(imax))
    print *,"tabulating dust emission integrals..."
    do j=1,nd !loop on bins
       do i=1,imax !loop on Tbb
          Tbb = (i-1) * (1d4-1d0) / (imax-1) + 1d0 !BB temperature
          dust_opt_Tbb(i) = Tbb !store Tbb values
          !store and compute integrals
          dust_opt_Em(j,i) = dustEm(krome_dust_asize(j),Tbb,&
               dust_opt_asize, dust_opt_nu, dust_opt_Qabs) 
       end do
    end do
    print *,"done!"
  end subroutine dustOptIntegral

  !**********************
  function getTdust(dust_opt_Tbb, dust_opt_Em, &
       dust_opt_asize, dust_opt_nu, dust_opt_Qabs,n)
    !Tdust evaluation from eqn.6 Schneider+2006 MNRAS_369
    use krome_commons
    use krome_constants
    integer::i,nd,j,j2,jr,jl,iter,nT
    real*8::asize,eAbs
    real*8::dustHl,dustHr,dustH,rhogr,sgnr,sgnl,dust2H
    real*8::ntot,n(:),xl,xr,x2,getTdust(ndust)
    real*8::dust_opt_Em(:,:),dust_opt_Tbb(:),dust_opt_asize(:)
    real*8::dust_opt_nu(:), dust_opt_Qabs(:,:), Tdust(ndust)
    logical::found
    nd = ndust
    ntot = sum(n(1:nmols)) !total density
    nT = size(dust_opt_Tbb) !number of temperature values
    !loop on bins
    do i=1,nd
       if(n(nmols+i)<1d-19) cycle !note that this limit depends on ATOL
       asize = krome_dust_asize(i) !mean bin size
       rhogr = n(nmols+i) * asize**3 * krome_grain_rho
       !absorbed erg/cm3/s (rhogr g/cm3)
       eAbs = dustAbs(asize, dust_opt_asize, dust_opt_nu, dust_opt_Qabs) * rhogr

       dustHl = eAbs - dust_opt_Em(i,1) *rhogr &
            + dustCool(krome_dust_asize2(i), n(nmols+i), n(idx_Tgas),&
            dust_opt_Tbb(1), ntot)
       dustHr = eAbs - dust_opt_Em(i,nT) *rhogr &
            + dustCool(krome_dust_asize2(i), n(nmols+i), n(idx_Tgas),&
            dust_opt_Tbb(nT), ntot)

       sgnr = sgn(dustHr) !stores right bound sign
       sgnl = sgn(dustHl) !stores left bound sign

       !checks for bound signs
       if(sgnr==sgnl) then
          print *,"ERROR: probably the opacity-temperature curve has no roots!"
          stop
       end if

       jl = 1
       jr = nT
       iter = 0
       !root-finding
       do
          iter = iter + 1
          j2 = .5d0 * (jr + jl)
          dust2H = eAbs - dust_opt_Em(i,j2) *rhogr &
               + dustCool(krome_dust_asize2(i), n(nmols+i), n(idx_Tgas),&
               dust_opt_Tbb(j2), ntot)
          if(sgn(dust2H)==sgnr) then
             jr = j2
          else
             jl = j2
          end if
          if(abs(jr-jl).le.1) exit
       end do
       dustHr = eAbs - dust_opt_Em(i,jr) *rhogr &
            + dustCool(krome_dust_asize2(i), n(nmols+i), n(idx_Tgas),&
            dust_opt_Tbb(jr), ntot)
       dustHl = eAbs - dust_opt_Em(i,jl) *rhogr &
            + dustCool(krome_dust_asize2(i), n(nmols+i), n(idx_Tgas),&
            dust_opt_Tbb(jl), ntot)
       Tdust(i) = (0.d0 - dustHr)/(dustHl-dustHr)&
            *(dust_opt_Tbb(jl)-dust_opt_Tbb(jr)) + dust_opt_Tbb(jr)
    end do
    getTdust(:) = Tdust(:)
  end function getTdust
  
  !*****************************
  function beta(n,asize,Tdust,nndust)
    !opacity term (Omukai+2000, 2005)
    use krome_commons
    use krome_constants
    use krome_subs
    real*8::n(:),ntot,beta,tau,lj,rhodust,Tdust
    real*8::nndust,asize,rhogas,Tgas,m(nspec)

    m(:) = get_mass() !get masses (g)
    Tgas = n(idx_Tgas) !get Tgas (K)
    rhodust = nndust * asize**3 * krome_grain_rho !dust mass density g/cm3
    rhogas = sum(n(1:nmols)*m(1:nmols)) !gas mass density g/cm3
    !jeans length
    lj = sqrt(15./4.*boltzmann_erg*Tgas/pi/gravity/1.22/rhogas/p_mass)
    ntot = sum(n(1:nmols)) !total density gas
    tau = kopa(Tdust) * rhogas * lj !opacity
    beta = min(1d0, tau**(-2)) !final opacity

  end function beta
  !****************************
  function kopa(Tbb)
    real*8::kopa,Tbb,x
    !returns black body integrated with 
    ! opacity as kernel in cm2/g (Dopcke+2012)
    x = Tbb
    
    if(x<2d2) then
       kopa = 0.0004d0*x*x
    elseif(x.ge.2d2 .and. x<1.5d3) then
       kopa = 16.d0
    else
       kopa = 2.07594d-12*x**(-12)
    end if

  end function kopa
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

 !*****************+
  function dustAbs(asize, dust_opt_asize, dust_opt_nu, dust_opt_Qabs)
    !find dust Absorption for a given size grain
    ! i.e. integrates J(nu)*Qabs(nu)*dnu with trapezoidal method
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_user_commons
    real*8::dustAbs,asize,int0,int1,dnu,kera,kerb
    real*8::dust_opt_asize(:),dust_opt_nu(:), dust_opt_Qabs(:,:)
    integer::i,ival

    !find the nearest (larger) size value to asize
    do i=1,size(dust_opt_asize)
       if(asize<dust_opt_asize(i)) then
          ival = i !store larger index (right)
          exit
       end if
    end do

    !intgrate for the two size values (left, right) close to asize
    ! then interpolate the two results
    int0 = 0.d0 !init left integral 
    int1 = 0.d0 !init right integral
    !computes integral by using trapezoidal method
!!$    do i=2,size(dust_opt_nu)
!!$       !flux is the kernel (J(nu)+Planck)
!!$       kera =  Jflux(dust_opt_nu(i-1)*planck_eV) * eV_to_erg + &
!!$            fplanck(dust_opt_nu(i-1), 2.73d0*(1.d0+krome_redshift)) 
!!$       kerb = Jflux(dust_opt_nu(i)*planck_eV) * eV_to_erg + &
!!$            fplanck(dust_opt_nu(i), 2.73d0*(1.d0+krome_redshift))
!!$       dnu = dust_opt_nu(ival) - dust_opt_nu(ival-1) !stepsize
!!$       !calculates integrals
!!$       int0 = int0 + 0.5d0 * (kera * dust_opt_Qabs(ival-1,i-1) &
!!$            + kerb * dust_opt_Qabs(ival-1,i)) * dnu
!!$       int1 = int1 + 0.5d0 * (kera * dust_opt_Qabs(ival,i-1) &
!!$            + kerb * dust_opt_Qabs(ival,i)) * dnu
!!$    end do

    !absorption: interpolates values depending on asize: erg/g/s
    dustAbs = (int1 - int0) * (asize - dust_opt_asize(ival-1)) &
         / (dust_opt_asize(ival) - dust_opt_asize(ival-1)) + int0
    
  end function dustAbs

  
  !*****************+
  function dustEm(asize, Tbb, dust_opt_asize, dust_opt_nu, dust_opt_Qabs)
    !find dust emittivity for a given size and black body temperature
    ! i.e. integrates B(nu)*Qabs(nu)*dnu with trapezoidal method
    use krome_commons
    real*8::dustEm,asize,int0,int1,dnu,Tbb,kera,kerb
    real*8::dust_opt_asize(:), dust_opt_nu(:), dust_opt_Qabs(:,:)
    integer::i,ival

    !find the nearest size value to asize
    do i=1,size(dust_opt_asize)
       if(asize<dust_opt_asize(i)) then
          ival = i !store index
          exit
       end if
    end do

    !intgrate for the two size values (left, right) close to asize
    ! then interpolate the two results
    int0 = 0.d0 !init left integral 
    int1 = 0.d0 !init right integral
    !computes integral by using trapezoidal method
    do i=2,size(dust_opt_nu)
       kera = fplanck(dust_opt_nu(i-1), Tbb) !Planck function is the kernel
       kerb = fplanck(dust_opt_nu(i), Tbb) !Planck function is the kernel
       dnu = dust_opt_nu(ival) - dust_opt_nu(ival-1) !stepsize
       !calculates integrals
       int0 = int0 + 0.5d0 * (kera * dust_opt_Qabs(ival-1,i-1) &
            + kerb * dust_opt_Qabs(ival-1,i)) * dnu
       int1 = int1 + 0.5d0 * (kera * dust_opt_Qabs(ival,i-1) &
            + kerb * dust_opt_Qabs(ival,i)) * dnu
    end do

    !emission: interpolates values depending on asize: erg/g/s
    dustEm = (int1 - int0) * (asize - dust_opt_asize(ival-1)) &
         / (dust_opt_asize(ival) - dust_opt_asize(ival-1)) + int0
    
    
  end function dustEm

  !***************
  subroutine init_Qabs(fname,opt_Qabs,opt_asize,opt_nu)
    use krome_commons
    integer::ios,icount,acount,nucount,i
    real*8::x(5),xold(size(x))
    real*8,allocatable::opt_Qabs(:,:),opt_asize(:),opt_nu(:)
    character(*)::fname
    print *,fname

    !skip reading file if already allocated
    if(allocated(opt_Qabs)) return
    
    !open file to count data lines
    open(71,file=fname,status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: problems with file gra.dat"
       stop
    end if

    icount = 0
    acount = 0
    nucount = 0
    x(:) = 0.d0
    do
       read(71,*,iostat=ios) x(:)
       if(ios<0) exit !eof
       if(ios.ne.0) then
          cycle !blank line
       end if
       if(xold(1).ne.x(1)) then
          acount = acount + 1
          if(nucount==0 .and. icount>1) nucount = icount
       end if
       icount = icount + 1
       xold(:) = x(:)
    end do
    close(71)
    
    !allocate Qabs array
    allocate(opt_Qabs(acount,nucount),opt_nu(nucount),&
         opt_asize(acount))

    print *,"acount:",acount
    print *,"nucount:",nucount
    print *,"totcount:",acount*nucount

    !load file data into Qabs array
    open(71,file=fname,status="old",iostat=ios)
    icount = 0
    acount = 0
    xold(:) = 0
    do
       read(71,*,iostat=ios) x(:)
       if(ios<0) exit !eof
       if(ios.ne.0) then
          cycle !blank line
       end if
       icount = icount + 1
       if(xold(1).ne.x(1)) then
          acount = acount + 1
          opt_asize(acount) = x(1)
       end if
       if(icount.le.size(opt_nu)) opt_nu(icount) = x(2)
       opt_Qabs(acount,mod(icount-1,nucount)+1) = x(3)
       xold(:) = x(:)
    end do
    close(71)

    !convert the Qabs to cm2/g
    do i=1,size(opt_Qabs,1)
       opt_Qabs(i,:) = opt_Qabs(i,:) / opt_asize(i) / krome_grain_rho 
    end do

    print *,"Dust optical properties loaded from: ", fname
        
  end subroutine init_Qabs

  !**********************
  function fplanck(nu,Tbb)
    !Planck law
    use krome_constants
    real*8::fplanck,nu,Tbb
    !nu: Hz
    !Tbb: K
    !fplanck: erg/cm2

    if(exp_planck*nu>1d2) then
       fplanck = 0.d0
       return
    end if
    fplanck = pre_planck * nu**3 / (exp(exp_planck*nu/Tbb) - 1d0) !erg/cm2
  end function fplanck


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
