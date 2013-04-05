program test_krome

  !***************************
  !This test simulates the evolution of a planet atmosphere
  !for different atmospheric layers
  !***************************

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_constants
  use IFPORT !remove when using GNU fortran compiler (gfortran)
  
  !nspec is shown during krome.py execution
  integer,parameter::nspec=57,nlayer=10
  integer::i,j
  real*8::x(nspec),Tgas,dt,t

  !loop over layers
  do i = 1,nlayer
     Tgas = 270.d0 !gas Temperature (K)
     x(:) = 0.d0 !initialize species to default
     !randomize initial species
     do j=1,nspec
        x(j) = rand()
     end do
     x(:) = x(:)/sum(x) * 1d10
     Pgas = 1.d1*(i*5./nlayer) !pressure layer dependent
     dt = seconds_per_year * 1d1 !initial time-step (s)
     t = 0.d0 !initial time (s)
     print *,"solving layer",i
     !loop over time to provide intermediate steps
     do 
        dt = dt * 2d0 !increase time-step (s)
        t = t + dt !increse time (s)
        call krome(x(:),Tgas,dt) !call KROME
        write(66,'(I5,99E12.3e3)') i,t/seconds_per_year, x(:) !write to file
        if(t>seconds_per_year*1d5) exit !break loop when overshoot 1e5 years
     end do
     write(66,*) !blank line at end of each loop
  end do
  print *,"done!"
  
  !bye bye
  print *,"that's all! have a nice day!"

end program test_krome
