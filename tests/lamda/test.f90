!###################################################

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use krome_user_commons ! include krome_invdvdz
  implicit none
  real*8::Tgas,x(krome_nmols)
  integer::i,imax
  
  
  call krome_init() !init krome (mandatory)
  
  x(:) = 1d0 !default abundances

  krome_invdvdz = 1d11
  imax = 30
  do i=1,imax
     Tgas = 1d1**((i-1)*(4d0-1d0)/(imax-1)+1d0)
     !x(:) = 1d1**((i-1)*(6d0-1d0)/(imax-1)+1d0)
     write(66,'(99E17.8e3)') Tgas, krome_coolingCI(x(:),Tgas), &
          krome_coolingCII(x(:),Tgas), &
          krome_coolingOI(x(:),Tgas), &
          krome_coolingCOI(x(:),Tgas), &
          krome_cooling13COI(x(:),Tgas)
  end do
    
  print *,"Test OK!"
  
end program test

