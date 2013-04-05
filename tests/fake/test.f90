program test_krome

  !********************************
  !This test uses fake species named FK1,FK2,... 
  !and a random-generated network of reactions.
  !The test is repeated imax times
  !with random-generated initial conditions.
  !********************************

  use krome_main
  use krome_subs
  use krome_user
  use krome_user_commons
  use IFPORT !REMOVE when using GNU compiler (gfortran)

  integer,parameter::nspec=3
  integer::i,imax,j
  real*8::x(nspec),Tgas,rho,dt,xold(nspec),t
  real*8::lmax,lmin

  lmax = log10(1.d0) !log of maximum initial abundances
  lmin = log10(1d-6) !log of minimum initial abundances
  rho = 10.d0 !gas density (g/cm3)
  Tgas = 0.826d2 !gas Temperature (K)
  imax = 30 !number of evolutions

  !loop over evolutions
  do i = 1,imax
     print '(a10,I5,a4,I5)',"Evolution",i,"of",imax
     !intialize species randomly
     do j=1,nspec
        x(j) = 1d1**(rand()*(lmax-lmin)+lmin)
     end do
     !normalize initial fractions (#)
     x(:) = x(:) / sum(x)
     dt = 1d3 !time-step (s)
     t = 0.d0 !initial time (s)
     write(66,'(999E12.3e3)') t,x(:)
     !loop over time to get intermediate integration steps
     do 
        xold(:) = x(:) !store fractions
        dt = dt * 1.5d0 !increase time-step
        call krome(x(:), rho, Tgas, dt) !call KROME
        t = t + dt !increase time
        write(66,'(999E12.3e3)') t,x(:)
        if(sum((xold(:)-x(:))**2)/nspec<1d-20) exit !check the steady-state
     end do
     write(77,'(999E12.3e3)') t,x(:)
     write(66,*)
     write(66,*)
  end do
  
  print *,"done!"
  print *,"Output written in fort.66"
  print *,"Default gnuplot: load 'plot.gps'"
  
  print *,"that's all! have a nice day!"

end program test_krome
