module krome_user
  implicit none

#KROME_header

#KROME_species

#KROME_common_alias

contains

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

  !*****************
  !alias for get_mass
  function krome_get_mass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec), krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
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
  !alias for get_index
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer::krome_get_index
    character*(*)::name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  function krome_get_rho(n)
    use krome_commons
    real*8::krome_get_rho,n(:)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho
  
  !*************************
  subroutine krome_scale_Z(n,Z)
    !scale metallicity according to Arnett(1996)
    use krome_commons
    real*8::n(:),Z

#KROME_scaleZ
    
  end subroutine krome_scale_Z

  !***********************
  function krome_get_electrons(n)
    use krome_commons
    real*8::n(:),ee,krome_get_electrons
    ee = sum(n(:) * krome_get_charges())
    krome_get_electrons = max(0.d0, ee)
  end function krome_get_electrons

  !*******************************
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


  !************************
  !print species informations
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
