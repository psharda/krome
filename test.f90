program smith
  use krome_main
  use krome_user
  implicit none
  real*8::x(krome_nmols),t,dt,Tgas,spy
  integer::j,jmax

  jmax = 30

  spy = krome_seconds_per_year
  call krome_init()

  do j=1,jmax

     Tgas = 1d1**((j-1)*(5d0-1d0)/(jmax-1)+1d0)

     x(:) = 0d0
     x(H2) = 1d0
     x(CO) = 1d0

     t = 0d0
     dt = spy
     call krome_thermo_OFF()
     do
        dt = dt * 1.1d0
        call krome(x(:),Tgas,dt)
        t = t + dt
        write(66,'(99E17.8e3)') t/spy,Tgas,x(:)
        if(t>1d8*spy) exit
     end do
     write(44,'(99E17.8e3)') ntot,Tgas,krome_cooling(x(:),Tgas)
  end do
  
end program smith
