module krome_user
  implicit none

#KROME_header

#KROME_species

contains

#IFKROME_use_cooling
  !******************
  subroutine plot_cooling(n)
    use KROME_cooling
    implicit none
    real*8::n(:)
    
    call plot_cool(n(:))

  end subroutine plot_cooling
#ENDIFKROME
end module krome_user
