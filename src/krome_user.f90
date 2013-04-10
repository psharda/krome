module krome_user
  implicit none

#KROME_header

#KROME_species

contains

#IFKROME_use_cooling
  !******************
  !alias of plot_cool
  subroutine plot_cooling(n)
    use KROME_cooling
    implicit none
    real*8::n(:)

    call plot_cool(n(:))

  end subroutine plot_cooling
#ENDIFKROME

  !*****************
  !alias for get_mass
  function krome_get_mass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::krome_get_mass(nspec)
    krome_get_mass = get_mass()
  end function krome_get_mass

  !*****************
  !alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    implicit none
    real*8::krome_get_charges(nspec)
    krome_get_charges = get_charges()
  end function krome_get_charges

  !*****************
  !alias for get_names
  function krome_get_names()
    use krome_subs
    use krome_commons
    implicit none
    character*16::krome_get_names(nspec)
    krome_get_names = get_names()
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
  
  
end module krome_user
