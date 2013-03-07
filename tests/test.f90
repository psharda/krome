program test_simple

  !this use statement is mandatory
  use krome_main
  !this use statement provide useful common variables
  ! like species index (e.g. KROME_idx_C for carbon
  ! or KROME_idx_H2 for molecular hydrogen H2)
  use krome_user
  implicit none
  
  !array of the species abundances the size is suggested by KROME
  ! at the end of the output
  real*8::x(5)
  !other variables
  real*8::temperature,density,time_step,time
  integer::i

  !this call initialize krome (mandatory)
  call krome_init()

  !set zero as default abundances (see also below)
  x(:) = 1.d-10

  !initialize chemical abundnaces not equal to zero (see above)
  ! as fraction of the total density
  x(krome_idx_CH) = 0.5d0
  x(krome_idx_H)  = 0.5d0
  
  temperature = 1d1 !gas temperature in K
  density = 1d-20 !gas density in g/cm3
  
  time = 0.d0 !initialize time in seconds
  time_step = 1.d-6 !timestep in seconds
  
  write(*,'(I4,999E11.3)') 0, time, x(:)
  !call krome several times to produce an evolving output
  do i=1,50
     !call the solver
     call krome(x(:), density, temperature, time_step)
     !increase the total time
     time = time + time_step
     time_step = time_step * 2.d0
     !write some output
     write(*,'(I4,999E11.3)') i, time, x(:)
  end do
end program test_simple
