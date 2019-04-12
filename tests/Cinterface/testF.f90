program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  integer:: nstep, unit
  integer, parameter:: ioin = 10
  integer, parameter:: nsp=krome_nmols !number of species (common)
  real*8:: Tgas, dt, x(nsp), spy, t, ntot, eVperg, fe60_xi, x_Fe60

  spy = krome_seconds_per_year !seconds per year

  call krome_init() !init krome (mandatory)

  ! ----==== Isotope decay off ====---- !

  x(:) = 1d-20 !default abundances
  ntot = 1.0d4 ! cm**-3
  x(krome_idx_Hj) = ntot !hydrogen ion initial abundance

  ! set iron abundance
  call krome_scale_Z(x(:),0.0d0)
  x(krome_idx_Fej) = x(krome_idx_Fe)
  x(krome_idx_Fe) = 0.0d0
  x(krome_idx_E) = krome_get_electrons(x(:))

  ! set 60Fe abundance
  call krome_set_user_tauh(1.5d6 * spy)
  x_Fe60 = 0.0d0
  print*, 'Fe-60 abundance =',x_Fe60
  x(krome_idx_60Fe) = x_Fe60 * x(krome_idx_Fej)
  ! set rate for ionisation by isotope decay
  fe60_xi = 0.0d0
  print*, 'Isotope decay off.'
  call krome_set_user_xi(fe60_xi)
  ! set heating from isotope decay
  call krome_set_user_Wergs(36.0d0 * krome_eV_to_erg)

  Tgas = 1.0d3 !gas temperature (K)

  open (newunit=unit, file='idoff.f.dat', status='unknown', form='formatted')

  t = 0.0d0
  dt = 1d-2 * spy ! yr -> s
  nstep = 0
  do
    dt = dt * 1.1d0
    t = t + dt
    !call KROME
    call krome(x(:), Tgas, dt)
    nstep = nstep + 1
    if (mod(nstep,10) .eq. 0) then
      write (*,'(a,i5.1)') "nstep = ", nstep
    endif
    write (unit,'(60(1pe15.8))') t/spy, Tgas, x(:)/ntot
    call krome_popcool_dump(t/spy, 69)
    if (t > 1.0d8 * spy) exit
  enddo
  close(unit)

  write (*,'(a,i5.1)') "Finished. Number of steps = ", nstep

  ! ----==== Isotope decay on ====---- !

  x(:) = 1d-20 !default abundances
  ntot = 1.0d4 ! cm**-3
  x(krome_idx_Hj) = ntot !hydrogen ion initial abundance

  ! set iron abundance
  call krome_scale_Z(x(:),0.0d0)
  x(krome_idx_Fej) = x(krome_idx_Fe)
  x(krome_idx_Fe) = 0.0d0
  x(krome_idx_E) = krome_get_electrons(x(:))
  ! set 60Fe abundance
  call krome_set_user_tauh(1.5d6 * spy)
  x_Fe60 =1.0d-6
  print*, 'Fe-60 abundance =',x_Fe60
  x(krome_idx_60Fe) = x_Fe60 * x(krome_idx_Fej)
  ! set rate for ionisation by isotope decay
  fe60_xi = 1.0d-10
  print*, 'Isotope decay on.'
  call krome_set_user_xi(fe60_xi)
  ! set heating from isotope decay
  call krome_set_user_Wergs(36.0d0 * krome_eV_to_erg)

  Tgas = 1.0d3 !gas temperature (K)

  open (newunit=unit, file='idon.f.dat', status='unknown', form='formatted')

  t = 0.0d0
  dt = 1d-2 * spy ! yr -> s
  nstep = 0
  do
    dt = dt * 1.1d0
    t = t + dt
    !call KROME
    call krome(x(:), Tgas, dt)
    nstep = nstep + 1
    if (mod(nstep,10) .eq. 0) then
      write (*,'(a,i5.1)') "nstep = ", nstep
    endif
    write (unit,'(60(1pe15.8))') t/spy, Tgas, x(:)/ntot
    call krome_popcool_dump(t/spy, 70)
    if (t > 1.0d8 * spy) exit
  enddo
  close(unit)

  write (*,'(a,i5.1)') "Finished. Number of steps = ", nstep

end program test
