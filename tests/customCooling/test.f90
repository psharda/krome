!***************
!this test shows how to use custom cooling functions
!check the network react_customCool for details
!***************
program test
  use krome_main
  use krome_user
  implicit none
  real*8::x(krome_nmols),Tgas
  integer::j,jmax

  jmax = 30

  call krome_init()

  write(44,'(a)') "#Tgas "//krome_get_cooling_names_header()
  !loop on Tgas points
  do j=1,jmax

     !compute Tgas
     Tgas = 1d1**((j-1)*(5d0-1d0)/(jmax-1)+1d0)

     x(:) = 0d0
     x(krome_idx_H2) = 1d0
     x(krome_idx_H2O) = 1d0
     x(krome_idx_H) = 1d0
     write(44,'(99E17.8e3)') Tgas,krome_get_cooling_array(x(:),Tgas)
  end do

  print *,"test DONE!"
  print *,"plot in gnuplot typing"
  print *, "gnuplot> load 'plot.gps'"

end program test
