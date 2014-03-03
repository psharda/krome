module krome_user
  implicit none

#KROME_header

#KROME_species

#KROME_common_alias

#KROME_constant_list

contains

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt)
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    implicit none
    real*8::x(:),Tgas,dt,dTgas,k(nrea),krome_gamma
    real*8::n(nspec),nH2dust

#IFKROME_use_thermo
    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
         * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas
#ENDIFKROME

  end subroutine krome_thermo


#IFKROME_use_cooling
  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n)
    use krome_cooling
    implicit none
    real*8::n(:)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nmols),Tgas,x(nspec)
    integer,optional::nfile_in
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
    if(present(nfile_in)) nfile = nfile_in
    call dump_cool(x(:),Tgas,nfile)
    
  end subroutine krome_dump_cooling

#ENDIFKROME

  !*************************
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8::x(:),xi(:)
    real*8::n(krome_nspec),ni(krome_nspec),krome_conserve(krome_nmols)
    
    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas)
    use krome_subs
    use krome_commons
    real*8::krome_get_gamma,x(nmols),n(nspec),Tgas
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies
  !alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    implicit none
    integer::krome_get_zatoms(nmols),zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)
    
  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from 
  ! number density and mass density
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8::krome_get_mu,x(:),rhogas,n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    rhogas = sum(x(:)*krome_get_mass())
    krome_get_mu = get_mu(n(:),rhogas)
  end function krome_get_mu

  !*****************
  !get an array of double containing the masses in g
  ! of the species
  !alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec), krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse 
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec), krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x)
    use krome_commons
    use krome_subs
    real*8::n(nspec),x(:),krome_get_Hnuclei
    n(:) = 0d0
    n(1:nmols) = x(:)
    
    krome_get_Hnuclei = get_Hnuclei(n(:))
    
  end function krome_get_Hnuclei

  !*****************
  !get an array containing the charges of the species
  !alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec),krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 containing the names
  ! of alla the species
  !alias for get_names
  function krome_get_names()
    use krome_subs
    use krome_commons
    implicit none
    character*16::krome_get_names(nmols),tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !*****************
  !get the index of the species with name name
  !alias for get_index
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer::krome_get_index
    character*(*)::name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n)
    use krome_commons
    real*8::krome_get_rho,n(:)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho
  
  !*************************
  subroutine krome_scale_Z(n,Z)
    !scale metallicity the metals contained in n(:) 
    ! to Z according to Arnett(1996)
    use krome_commons
    real*8::n(:),Z

#KROME_scaleZ
    
  end subroutine krome_scale_Z

  !***********************
  function krome_get_electrons(n)
    !get the total number of electrons from
    ! the number denisities of all the species
    use krome_commons
    real*8::n(:),ee,krome_get_electrons,x(size(n))
    x(:) = n(:)
#KROME_zero_electrons
    ee = sum(x(:) * krome_get_charges())
    krome_get_electrons = max(0.d0, ee)
  end function krome_get_electrons


  !**********************
  !print the nbest fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest)
    use krome_subs
    use krome_commons
    implicit none
    real*8::x(nmols),xin(nmols),n(nspec),Tgas
    integer::nbest
    x(:) = xin(:)
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)
  end subroutine krome_print_best_flux

  !*******************************
  !get the fluxes of all the reactions in 1/cm3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8::krome_get_flux(nrea),n(nmols),Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*********************
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    implicit none
    real*8::krome_get_qeff(nrea)
    
    krome_get_qeff(:) = get_qeff()
    
  end function krome_get_qeff

#IFKROME_useStars

  !**************************
  !alias for stars_coe in krome_star
  function krome_stars_coe(x,rho,Tgas,y_in,zz_in)
    use krome_commons
    use krome_stars
    use krome_subs
    implicit none
    real*8::rho,Tgas,krome_stars_coe(nrea)
    real*8::n(nspec),x(:)
    integer::ny
    real*8,allocatable::y(:)
    integer,allocatable::zz(:)
    real*8,optional::y_in(:)
    integer,optional::zz_in(:)
    
    !check if extened abundances and zatom array are present
    if(present(y_in)) then
       ny = size(y_in)
       allocate(y(ny), zz(ny))
       y(:) = y_in(:)
       zz(:) = zz_in(:)
    else
       allocate(y(nmols),zz(nmols))
       zz(:) = get_zatoms()
       y(:) = x(:)
    end if

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    
    krome_stars_coe(:) = stars_coe(n(:),rho,Tgas,y(:),zz(:))
    deallocate(y,zz)
  end function krome_stars_coe

  !********************************
  !alias for stars_energy
  function krome_stars_energy(x,rho,Tgas,k)
    use krome_commons
    use krome_stars
    implicit none
    real*8::x(:),rho,Tgas,krome_stars_energy(nrea)
    real*8::n(nspec),k(nrea)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_stars_energy(:) = stars_energy(n(:),rho,Tgas,k(:))

  end function krome_stars_energy
    
#ENDIFKROME

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile)
    use krome_commons
    real*8::n(:),Tgas,flux(nrea)
    integer::nfile,i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
       write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas)
    use krome_commons
    use krome_subs
    integer::i,charges(nspec)
    real*8::x(:),Tgas,masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
       print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info
  
end module krome_user
