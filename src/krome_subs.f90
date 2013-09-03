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
    real*8::Tgas2,Tgas3,Tgas4,T0,T02,T03,T04,T0inv
    integer::i
    !Tgas is in K
    
#KROME_Tshortcuts

    k(:) = 1.d-40 !inizialize coefficients

#KROME_krates

    !set a minumum value for the rates
    !do i=1,nrea
    !   k(i) = max(k(i),1d-40) 
    !end do

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
