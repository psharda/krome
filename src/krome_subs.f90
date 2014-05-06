module krome_subs
contains

#KROME_header

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_user_commons
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),t
#KROME_shortcut_variables
    real*8::small,nmax,fsh
    integer::i
#KROME_initcoevars
    !Tgas is in K
    Tgas = max(n(idx_Tgas), 2.73d0)
    T = Tgas

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = #KROME_small

#KROME_Tshortcuts

#KROME_coevars

#IFKROME_useShieldingDB96
    !compute shielding from Draine+Bertoldi 1996
    fsh = calc_H2shieldDB96(n(:), Tgas)
#ENDIFKROME

#IFKROME_useShieldingWG11
    !compute shielding from Wolcott+Greene 2011
    fsh =  calc_H2shieldWG11(n(:), Tgas)
#ENDIFKROME

    k(:) = small !inizialize coefficients

#KROME_krates

    coe(:) = k(:) !set coefficients to return variable

  end function coe

#KROME_metallicity_functions  

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
    b5 = b*10d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5) + &
         0.035d0/sqrt(1d0+x) * &
         exp(-8.5d-4*sqrt(1+x))

  end function fselfH2
  
  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas_in,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas,Tgas_in
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    !avoid low-temperature problems
    Tgas = max(Tgas_in,1d1)
    
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
  function zf(Tgas_in,brot)
    real*8::Tgas,zf,brot,z,ibTgas,Tgas_in
    integer::j,jmax
    jmax = 10 !number of levels
    
    !avoid low-temperature problems
    Tgas = max(Tgas_in,1d1)

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
  function gamma_rotop(Tgas,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

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
  function gamma_rot(Tgas,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

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

    !get_mu = 1.22d0

  end function get_mu

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

#KROME_masses

  end function get_mass

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
    rhogas = sum(n(1:nmols)*m(1:nmols))
    mu = 1.22d0
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
         /p_mass/gravity/mu)

  end function get_jeans_length

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
       print '(I4,a1,a50,E17.8)',i," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

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

end module krome_subs
