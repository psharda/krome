module krome_subs
contains

#KROME_header

  !************************
  !compute reaction rates
  function coe(n)
    use krome_commons
    use krome_user_commons
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,t,t3,n(nspec)
    real*8::logT,lnT,Te,lnTe,T32,invT,invTe,sqrTgas,invsqrT32,sqrT32
    integer::i
    !Tgas is in K
    Tgas = max(2.73d0, n(idx_Tgas))
    t = Tgas !alias for Tgas (K)
    logT = log10(Tgas) !log10 of Tgas (#)
    lnT = log(Tgas) !ln of Tgas (#)
    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    T32 = Tgas/3.d2 !Tgas/(300 K) (#)
    t3 = T32 !alias for T32 (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    invTe = 1.d0/Te !inverse of T (1/eV)
    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)
    invsqrT32 = 1.d0/sqrt(T32)
    sqrT32 = sqrt(T32)

    k(:) = 0.d0 !inizialize coefficients

#KROME_krates

    !set a minumum value for the rates
    do i=1,nrea
       k(i) = max(k(i),1d-30) 
    end do

    coe(:) = k(:)!set coefficients to return variable
  end function coe

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

#KROME_masses

  end function get_mass

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

#KROME_names

  end function get_names

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


  !************************
  !print species informations
  subroutine get_infos(x, Tgas)
    use krome_commons
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
  end subroutine get_infos

  !*****************************
  subroutine load_arrays()
    use krome_commons

#KROME_implicit_arrays

  end subroutine load_arrays

end module krome_subs
