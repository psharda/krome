!################################################################
! This test is designed to check if CO, OH, H2O, and HCN cooling
! produce reasonable temperature evolutions. WARNING: do not use
! the network for production
!###############################################################
program test_krome
  use krome_main
  use krome_user
  use krome_user_commons
  implicit none
  real*8::ntot,Tgas,x(krome_nmols),dt,t,tend
  integer::unit

  ! init krome (load network tables)
  call krome_init()

  ! initial conditions
  ntot = 1d6  ! total density, 1/cm3
  Tgas = 3d2  ! temperature, K
  tend = 1d8 * krome_seconds_per_year ! integration time, s

  !species initialization, 1/cm3
  x(:) = 1d-40
  x(KROME_idx_H)  = ntot
  x(KROME_idx_H2) = 1d-6 * ntot
  x(KROME_idx_CN) = 7d-5 * ntot
  x(KROME_idx_HCO) = 1d-4 * ntot
  x(KROME_idx_O) = 1d-4 * ntot

  ! set cosmic rays
  call krome_set_user_crate(1.3d-17)

  ! open file to write output
  open(newunit=unit, file="data.out", status="replace")

  t = 0d0
  dt = krome_seconds_per_year

  ! pseudo time-dependent loop
  do
     dt = dt * 1.1
     ! solve the chemistry
     call krome(x(:), Tgas, dt)
     t = t + dt
     write(unit,'(99E17.8e3)') t / krome_seconds_per_year, Tgas, x(:) / ntot
     if(t>=tend) exit
  end do

  close(unit)

  print *, "To plot output use gnuplot as"
  print *, "gnuplot> load 'plot.gps'"
  print *, "That's all! have a nice day!"

end program test_krome
