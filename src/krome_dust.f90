module krome_dust

#IFKROME_useDust
contains
  subroutine krome_init_dust(xdust,adust,ntot,alow_arg,aup_arg,phi_arg)
    !krome_init_dust: initialize the dust ditribution (xdust)
    ! and the dust bin mean sizes (adust). Arguments are
    ! ntot(number_of_dust_types) the total abundance per 
    ! each bin size, alow_arg the size of the smallest
    ! bin size, aup_arg the largest, phi_arg the exponent
    ! of the MRN power law.
    use krome_commons
    use krome_subs
    implicit none
    real*8,optional::alow_arg,aup_arg,phi_arg
    real*8::iphi1,c,phi1,abin(ndust+1),mass(nspec),xdust(ndust),myc
    real*8::alow,aup,phi,ntot(ndustTypes),adust(ndust),Tbb,myx,a0,a1
    integer::i,j,ilow,iup,imax,nd

#KROME_dustPartnerIndex

    !default values
    alow = 5d-7 !lower size (cm)
    aup = 2.5d-5 !upper size (cm)
    phi = -3.5d0 !MNR distribution exponent (with its sign)
    phi1 = phi + 1.d0
    iphi1 = 1.d0/phi1
    if(present(alow_arg)) alow = alow_arg
    if(present(aup_arg)) aup = aup_arg
    if(present(phi_arg)) phi = phi_arg
    mass(:) = get_mass()
    nd = ndust/ndustTypes

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
       krome_dust_asize(ilow:iup) = adust(:) !store mean size
       krome_dust_asize2(ilow:iup) = adust(:)**2 !store mean size squared
       xdust(ilow:iup) = ntot(j) / nd !amount of dust per bin

       !evaluate dust-parnter ratio (e.g. 1dust=1e2 C atoms)
       krome_dust_partner_ratio(ilow:iup) = adust(ilow:iup)**3 &
            * 2.3d0  / mass(krome_dust_partner_idx(j))
       krome_dust_partner_ratio_inv(ilow:iup) = 1.d0 &
            / krome_dust_partner_ratio(ilow:iup)
    end do
    

    !compute the mass of the dust partner
    do j=1,ndustTypes
       krome_dust_partner_mass(j) =  mass(krome_dust_partner_idx(j))
    end do

    krome_dust_T(:) = 5.d0 !defualt dust temperature

    !init optical properties
#KROME_init_Qabs

    !init integral Qabs(nu)*B(nu)*dnu
#KROME_opt_integral
    
    !set dust temperature
#KROME_getTdust

    print *,"Dust initialized!"

  end subroutine krome_init_dust


  !***********************
  subroutine dustOptIntegral(dust_opt_Em, dust_opt_Tbb, dust_opt_asize, &
       dust_opt_nu, dust_opt_Qabs)
    use krome_commons
    integer::i,imax,j,nd
    real*8,allocatable::dust_opt_Em(:,:), dust_opt_Tbb(:)
    real*8::dust_opt_asize(:), dust_opt_nu(:), dust_opt_Qabs(:,:),Tbb

    nd = size(krome_dust_asize)
    !computes and stores Qabs(nu)*B(nu)*dnu integral for each dust bin
    imax = int(1e3) !number of Tbb values
    allocate(dust_opt_Em(nd,imax), dust_opt_Tbb(imax))
    do j=1,nd !loop on bins
       do i=1,imax !loop on Tbb
          Tbb = (i-1) * (1d2-1d0) / (imax-1) + 1d0 !BB temperature
          dust_opt_Tbb(i) = Tbb !store Tbb values
          dust_opt_Em(j,i) = dustEm(krome_dust_asize(j),Tbb,&
               dust_opt_asize, dust_opt_nu, dust_opt_Qabs) !store and compute integrals
       end do
    end do
  end subroutine dustOptIntegral

  !**********************
  subroutine getTdust(Tdust, dust_opt_Tbb, dust_opt_Em, &
       dust_opt_asize, dust_opt_nu, dust_opt_Qabs)
    !Tdust evaluation from eqn.6 Schneider+2006 MNRAS_369
    use krome_commons
    use krome_constants
    integer::i,nd,j,nT
    real*8::asize,eAbs,T0,T1,mysgn,de0,de1,Tdust(:)
    real*8::dust_opt_Em(:,:),dust_opt_Tbb(:),dust_opt_asize(:)
    real*8::dust_opt_nu(:), dust_opt_Qabs(:,:)
    logical::found
    nd = ndust/ndustTypes !dust bins
    nT = size(dust_opt_Tbb) !number of temperature values
    !loop on bins
    do i=1,nd
       asize = krome_dust_asize(i) !mean bin size
       eAbs = dustAbs(asize, dust_opt_asize, dust_opt_nu, dust_opt_Qabs) !absorbed erg/cm2/s
       found = .false. !flag interpolation
       mysgn = sgn(eAbs - dust_opt_Em(i,1)) !first sign
       !loop on temperatures (discrete root-finding) 
       do j=2,nT
          de0 = dust_opt_Em(i,j-1) - eAbs !left function
          de1 = dust_opt_Em(i,j) - eAbs !right function
          !get change of sign
          if(mysgn * de1 > 0.d0) then
             !interpolate Tdust value
             Tdust(i) = (0.d0 - de0) / (de1 - de0) * (dust_opt_Tbb(j) - dust_opt_Tbb(j-1)) &
                  + dust_opt_Tbb(j-1)
             found = .true. !found flag
             exit
          end if
       end do
       
       !check found flag
       if(not(found)) then
          print *,"ERROR: dustEm not found!!!"
          stop
       end if
    end do
  end subroutine getTdust

  !******************
  function sgn(arg)
    real*8::sgn,arg
    sgn = 1.d0
    if(arg>=0.d0) return
    sgn = -1.d0
  end function sgn


  !*********************
  function krome_dust_grow(ndust,natom,Tgas,Tdust,vgas,adust)
    !krome_dust_grow: compute dust formation in cm3/s 
    ! (Grassi2012, eqn.25)
    use krome_constants
    implicit none
    real*8::krome_dust_grow,ndust,natom,Tgas,Tdust,vgas,adust,seed
    
    seed = 1d-12 !1/cm3
    krome_dust_grow = pi * adust**2 * (max(ndust,0.d0) + seed) * max(natom,0.d0) &
         * krome_dust_stick(Tgas,Tdust) * vgas
        
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
  function krome_dust_sput_DS79(Tgas,adust,natom,ndust)
    !krome_dust_sput_D97: sputtering following (Draine et al. 79)
    real*8::krome_dust_sput_DS79,Tgas,adust,natom,T6,ndust
    
    if(Tgas<1d4) then
       krome_dust_sput_DS79 = 0.d0
       return
    end if
    T6 = 1d6/Tgas
    krome_dust_sput_DS79 = 3.d-17 / adust * natom / (1.d0 + T6**3) * ndust !1/s

    if(krome_dust_sput_DS79>1.d0) then
       print *,krome_dust_sput_DS79,adust,natom,Tgas,T6,ndust
       stop
    end if
    !print *,natom,Tgas,T6,adust,krome_dust_sput_DS79
  end function krome_dust_sput_DS79

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
    do i=2,size(dust_opt_nu)
       kera =  Jflux(dust_opt_nu(i-1)*planck_eV) * eV_to_erg + &
            fplanck(dust_opt_nu(i-1), 2.73d0*(1.d0+redshift)) !flux is the kernel
       kerb = Jflux(dust_opt_nu(i)*planck_eV) * eV_to_erg + &
            fplanck(dust_opt_nu(i), 2.73d0*(1.d0+redshift)) !flux is the kernel
       dnu = dust_opt_nu(ival) - dust_opt_nu(ival-1) !stepsize
       !calculates integrals
       int0 = int0 + 0.5d0 * (kera * dust_opt_Qabs(ival-1,i-1) &
            + kerb * dust_opt_Qabs(ival-1,i)) * dnu
       int1 = int1 + 0.5d0 * (kera * dust_opt_Qabs(ival,i-1) &
            + kerb * dust_opt_Qabs(ival,i)) * dnu
    end do

    !absorption: interpolates values depending on asize: erg/cm2/s
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

    !emission: interpolates values depending on asize: erg/cm2/s
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
  function krome_H2_dust(ndust,Tdust,n,H2_eps_f,myvgas)
    !H2 formed on dust (1/cm3/s)
    use krome_constants
    use krome_commons
    real*8::H2_dust, krome_H2_dust,Tgas,Tdust(:)
    real*8::myvgas,H2_eps,ndust(:),n(:),H2_eps_f
    integer::i

    Tgas = n(idx_Tgas)
    
    H2_dust = 0.d0 
    do i = 1,size(Tdust)
       H2_eps = H2_eps_f(Tgas, Tdust(i))
       H2_dust = H2_dust + 0.5d0 * n(idx_H) * myvgas * ndust(i) * krome_dust_asize(i)**2 &
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
    apc = 1.7d0 * 1d-10 !m (must be in m even if in CS2009 it's in \AA!!, Cazaux2012 private comm.)
    func = 2.d0 * exp(-(Ep-Es)/(Ep+myTgas)) / (1.d0+sqrt((Ec-Es)/(Ep-Es)))**2
    H2_eps_Si = 1.d0/(1.+16.*myTdust/(Ec-Es)) * exp(-Ep/myTdust)&
          * exp(-4d9*apc*sqrt(Ep-Es)) + func

    if(H2_eps_Si>1.d0 .or. H2_eps_Si<0.d0) then
       print *,"problem on H2_eps_Si"
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
    H2_eps_C = (1.d0-TH) / (1.d0+0.25*(1.d0+sqrt((Ec-Es)/(Ep-Es)))**2 * exp(-Es/myTdust))
    if(H2_eps_C>1.d0 .or. H2_eps_C<0.d0) then
       print *,"problem on H2_eps_C"
       stop
    end if
  end function H2_eps_C

  !***************************
  function stick(myTgas,myTdust)
    !sticking coefficient per la formazione di H2
    real*8::stick,myTgas,myTdust
    stick = 1.d0 / (1.d0 + 0.4*sqrt((myTgas+myTdust)/1d2) + 0.2*myTgas/1d2 + 0.08*(myTgas/1d2)**2)
    if(stick>1.d0 .or. stick<0.d0) then
       print *,"problem on stick coefficient"
       stop
    end if
  end function stick



  
#ENDIFKROME
end module krome_dust
