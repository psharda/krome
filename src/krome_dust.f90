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
    real*8::iphi1,c,phi1,abin(ndust+1),mass(nspec),xdust(ndust)
    real*8::alow,aup,phi,ntot(ndustTypes),adust(ndust),nd
    integer::i,j,ilow,iup

#KROME_dustPartnerIndex

    !default values
    alow = 5d-8 !lower size (cm)
    aup = 2.5d-5 !upper size (cm)
    phi = -3.5d0 !MNR distribution exponent (with its sign)
    if(present(alow_arg)) alow = alow_arg
    if(present(aup_arg)) aup = alow_arg
    if(present(phi_arg)) phi = phi_arg
    mass(:) = get_mass()
    nd = ndust/ndustTypes
    !loop on dust types
    do j=1,ndustTypes
       ilow = nd * (j - 1) + 1 !lower index
       iup = nd * j !upper index
       phi1 = phi + 1.d0 !phi+1
       iphi1 = 1.d0/phi1 !1/phi
       
       !estimate normalization constant
       c = (phi1)/(aup**(phi1)-alow**(phi1))
       abin(1) = alow !set lower limit for the first bin
       !evaluates bin limits
       do i=2,nd+1
          abin(i) = (phi1/nd/c + abin(i-1)**phi1)**iphi1
       end do
       !evaluate means size of dust in each bin
       do i=1,nd
          adust(i+ilow-1) = abs(abin(i)-abin(i+1))*0.5d0
       end do
       !compute and normalize bin distribution
       xdust(ilow:iup) = c*adust(ilow:iup)**phi
       xdust(ilow:iup) = xdust(ilow:iup)/sum(xdust(ilow:iup)) * ntot(j)
       krome_dust_asize(ilow:iup) = adust(:)
       !evaluate dust-parnter ratio (e.g. 1dust=1e2 C atoms)
       krome_dust_partner_ratio(ilow:iup) = adust(ilow:iup)**3 &
            * 2.3d0  / mass(krome_dust_partner_idx(j))
       krome_dust_partner_ratio_inv(ilow:iup) = 1.d0 &
            / krome_dust_partner_ratio(ilow:iup)
    end do

    krome_dust_T(:) = 25.d0 !K

    !compute the mass of the dust partner
    do j=1,ndustTypes
       krome_dust_partner_mass(j) =  mass(krome_dust_partner_idx(j))
    end do
    print *,"Dust initialized!"

  end subroutine krome_init_dust

  !*********************
  function krome_dust_grow(ndust,natom,Tgas,Tdust,vgas,adust)
    !krome_dust_grow: compute dust formation in cm3/s 
    ! (Grassi2012, eqn.25)
    use krome_constants
    implicit none
    real*8::krome_dust_grow,ndust,natom,Tgas,Tdust,vgas,adust
    
    krome_dust_grow = pi * adust**2 * max(ndust,0.d0) * max(natom,0.d0) &
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
  function yield(energy_erg,mdust,matom,zidust,ziatom,kfree,U0_eV)
    !yield: dust sputtering from Nozawa2006.
    ! energy: impact energy (erg)
    ! matom: mass projectile atom (amu) e.g. H=1
    ! mdust: mass target dust (amu) e.g. C-dust=6
    ! zdust: atomic number dust, e.g. C-dust=12
    ! zatom: atomic number projectile, e.g. H=1
    ! kfree: free parameter see Nozawa2006
    ! U0_eV: see Nozawa (eV)
    implicit none
    real*8::yield,U0,kfree,mui,eth,gi,energy_erg
    real*8::mdust,matom,U0_eV,energy,zdust,zatom
    integer::zidust,ziatom
    mui = mdust/matom !mass ratio
    zdust = zidust
    zatom = ziatom
    U0 = U0_eV * 1.60217646d-12 !eV->erg
    energy = energy_erg !/ 1.60217646d-12 !erg->eV

    if(matom.le.0.3*mdust) then
       gi = 4.d0*matom*mdust/(matom+mdust)**2 !max transfer energy head-on
       eth = U0/gi*(1.d0-gi) !threshold energy
    else
       eth = 8.d0*U0*(matom/mdust)**(1./3.)
    end if
    if(eth/energy>1.) then
       print *,eth,energy,U0
       stop
    end if
    yield = 4.2d14 * nuclear_stop(energy,mdust,matom,zdust,zatom) / U0 *&
         alpha(mui)/(kfree*mui+1.d0) *&
         (1.d0-(eth/energy)**(2./3.)) *&
         (1.d0-eth/energy)**2

  end function yield

  !******************
  function alpha(mui)
    real*8::alpha,mui
    if(mui.le.0.5d0) then
       alpha = 0.2d0
    elseif(mui>1.d0) then
       alpha = 0.3d0*(mui-0.6)**(2./3.)
    else
       alpha = 0.1d0/mui + 0.25d0*(mui-0.5d0)**2
    end if
  end function alpha

  !******************
  function nuclear_stop(energy,mdust,matom,zdust,zatom)
    real*8::nuclear_stop,energy,eps,asc,pi,a0,ecgs
    real*8::mdust,matom,zdust,zatom
    a0 = 0.529d-8 !bohr raidus cm
    pi = 3.14159d0 
    ecgs = 4.80320425d-10 !electron charge (statC)
    asc = 0.855*a0/sqrt(zdust**(2./3.) + zatom**(2./3.)) !screening length
    eps = mdust/(matom+mdust)*asc/zdust/zatom/ecgs**2*energy !reduced energy
    nuclear_stop = 4.d0*pi*asc*zatom*zdust*ecgs**2*matom/(matom+mdust)*&
         columb_screen(eps) !nuclear stopping cross-section
  end function nuclear_stop

  !*******************
  function columb_screen(eps)
    real*8::columb_screen,eps
    columb_screen = (3.441d0*sqrt(eps)*log(eps+2.718d0)) / &
         (1.d0 + 6.35d0*sqrt(eps)+eps*(-1.708d0+6.882d0*sqrt(eps)))
  end function columb_screen

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


  !******************
  function krome_dust_sputtering(vgas,zatom,zdust,natom,ndust,adust,Tgas)
    !krome_dust_sputtering: compute dust sputtering in cm3/s
    ! from Nozawa2006
    use krome_constants
    real*8::krome_dust_sputtering,vgas,adust,Tgas,y,energy_1d4,y_1d4
    real*8::matom,mdust,kfree,U0_eV,ndust,natom(:),sput,energy
    integer::i,zdust,zatom(:)
        
    !sputtering is efficient only for Tgas>1d5
    if(Tgas<1d6) then
       krome_dust_sputtering = 0.d0
       return
    end if

    energy = Tgas * boltzmann_erg !energy at Tgas K
    energy_1d4 = 1d4 * boltzmann_erg !energy at 1d4 K

    if(zdust==12) then
       mdust = 12.d0
       U0_eV = 4.d0
       kfree = 0.61d0
    elseif(zdust==28) then
       mdust = 28.d0
       U0_eV = 4.66d0
       kfree = 0.43d0
    elseif(zdust==56) then
       mdust = 56.d0
       U0_eV = 4.31d0
       kfree = 0.23d0
    else
       print *,"ERROR: dust atomic number unknown!"
       print *,"Z number:",zdust
       stop
    end if

    sput = 0.d0
    do i=1,size(zatom)
       if(ndust.le.0.d0) cycle
       if(natom(i).le.0.d0) cycle
       matom = zatom(i)  !WARNING: rough approx. for atom mass
       y = yield(energy, mdust, matom, zdust, zatom(i), kfree, U0_eV)
       y_1d4 = yield(energy_1d4, mdust, matom, zdust, zatom(i), kfree, U0_eV)
       if(y<0.d0 .or. y>1.d0) then
          print *,"yield @ Tgas <0 or >1"
          print *,y
          stop
       end if
       if(y_1d4<0.d0 .or. y_1d4>1.d0) then
          print *,"yield @ 1d4 K <0 or >1"
          print *,y_1d4
          stop
       end if
       sput = sput + 1d-10 * ndust / adust * y / y_1d4
    end do
    print *,sput
    krome_dust_sputtering = sput
  end function krome_dust_sputtering

#ENDIFKROME
end module krome_dust
