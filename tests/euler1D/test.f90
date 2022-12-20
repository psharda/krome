! THIS MODULE EVOLVES 1D HYDRODYNAMICAL EQUATIONS USING HLL+RK2
! CLASSIC SOD SHOCK TEST, 0.125 density factor, 0.1 pressure factor
! NORMALIZED TO CGS UNITS
module hydro
  real*8,parameter::pc2cm=3.08567758128e18  ! cm / pc
  real*8,parameter::kboltzmann=1.380649e-16  ! erg/K
  real*8,parameter::pmass=1.67262192e-24  ! g

  integer,parameter::ncell=1024  ! number of grid cells
  integer,parameter::nall=ncell+2  ! number of grid cells including ghosts
  integer,parameter::idx_rho=1  ! index of density in hydro variables array, g/cm3
  integer,parameter::idx_rvx=2  ! index of momentum in hydro variables array, g/cm2/s
  integer,parameter::idx_energy=3  ! index of energy density in hydro variables array, erg/cm3
  integer,parameter::nhydro=3  ! number of hydro variables
  integer,parameter::nscalar=12  ! number of advected scalars (chemical species = krome_nmols)
  integer,parameter::nvars=nhydro + nscalar  ! total number of hydro variables
  real*8,parameter::gamma=5./3.  ! adiabatic index
  real*8,parameter::mu=1.25  ! mean molecular weight (fixed during evolution)
  real*8,parameter::dx=pc2cm / ncell  ! grid spacing, cm
contains

  ! ************************
  ! RHS of hydro equations
  function rhs(u) result(dy)
    implicit none
    real*8,intent(in)::u(ncell, nvars)
    real*8::ug(nall, nvars), p(nall)
    real*8::fg(nall, nvars), vx(nall), cs(nall)
    real*8::dy(ncell, nvars), lp(nall), lm(nall), ap, am
    real*8::fp(ncell, nvars), fm(ncell, nvars)
    integer::i

    ! set ghost cells
    ug(2:nall-1, :) = u(:, :)
    ug(1, :) = u(1, :)
    ug(nall, :) = u(ncell, :)

    ! compute velocity and pressure for later use
    vx = ug(:, idx_rvx) / ug(:, idx_rho)
    p = (ug(:, idx_energy) - ug(:, idx_rho) * vx**2 / 2.) * (gamma - 1d0)

    ! fluxes HD from conservation laws
    fg(:, idx_rho) = ug(:, idx_rvx)
    fg(:, idx_rvx) = ug(:, idx_rvx) * vx + p
    fg(:, idx_energy) = (ug(:, idx_energy) + p) * vx

    ! chemical species advection, g/cm3
    do i=nhydro+1,nvars
      fg(:, i) = vx * ug(:, i)
    end do

    ! speed of sound, cm/s
    cs = sqrt(gamma * p / ug(:, idx_rho))

    ! lambdas
    lp = vx + cs
    lm = vx - cs

    ! compute left and right fluxes
    do i=2,nall-1
      ap = max(0d0, max(lp(i), lp(i+1)))
      am = max(0d0, max(-lm(i), -lm(i+1)))
      fp(i-1, :) = ap * fg(i, :) + am * fg(i+1, :) - ap * am * (ug(i+1, :) - ug(i, :))
      fp(i-1, :) = fp(i-1, :) / (ap + am)

      ap = max(0d0, max(lp(i), lp(i-1)))
      am = max(0d0, max(-lm(i), -lm(i-1)))

      fm(i-1, :) = ap * fg(i-1, :) + am * fg(i, :) - ap * am * (ug(i, :) - ug(i-1, :))
      fm(i-1, :) = fm(i-1, :) / (ap + am)
    end do

    ! spatial gradient (RHS)
    dy = - (fp - fm) / dx

  end function rhs

  ! ********************
  ! use RK2 scheme to adavence the solution
  subroutine solve(rho, vx, tgas, xall, tend)
    use krome_main
    implicit none
    real*8,intent(inout)::rho(ncell), vx(ncell), tgas(ncell), xall(ncell, nscalar)
    real*8,intent(in)::tend
    real*8::u(ncell, nvars), p(ncell), cs(ncell), y0(ncell, nvars)
    real*8::k2(ncell, nvars), dt, t
    real*8::u1(ncell, nvars), cmax, vmax, dtmax, tgas1, p1, x(nscalar)
    integer::i
    cmax = 1d-1  ! Courant factor
    dtmax = dx / 1d6  ! largest time-step allowed, s

    ! copy initial conditions for all cells from subroutine arguments
    p = rho * kboltzmann * tgas / mu / pmass  ! pressure, erg/cm3

    ! speed of sound, cm/s
    cs = sqrt(gamma * p / rho)

    u(:, idx_rho) = rho  ! density, g/cm3
    u(:, idx_rvx) = rho * vx  ! momentum, g/cm2/s
    u(:, idx_energy) = p / (gamma - 1d0) + rho * vx**2 / 2.  ! energy density, erg/cm3
    u(:, nhydro+1:nvars) = xall(:, :)  ! chemical species, g/cm3

    y0 = u
    t = 0d0
    do
      ! compute time-step
      vmax = maxval(abs(y0(:, idx_rvx) / y0(:, idx_rho))+cs) + 1d-40
      vmax = max(vmax,maxval(abs(y0(:, idx_rvx) / y0(:, idx_rho))-cs) + 1d-40)
      dt = min(dx / vmax * cmax, dtmax)
      dt = min(dt, tend - t)

      ! advance RK2
      u1 = y0 + 0.5 * dt * rhs(y0)
      k2 = dt * rhs(u1)
      y0 = y0 + k2

      ! loop on cells to do call KROME
      do i = 1, ncell
         x = y0(i, nhydro+1:nvars) / y0(i, idx_rho)  ! mass fractions
         p1 = (y0(i, idx_energy) - y0(i, idx_rvx)**2 / y0(i, idx_rho) / 2.) * (gamma - 1d0)  ! pressure, erg/cm3
         tgas1 = p1 / y0(i, idx_rho) / kboltzmann * mu * pmass  ! tgas, K
         call krome(x, y0(i, idx_rho), tgas1, dt)  ! CALL KROME HERE
         y0(i, nhydro+1:nvars) = x * y0(i, idx_rho)  ! update abundances, g/cm3
         p1 = tgas1 / mu / pmass * kboltzmann * y0(i, idx_rho)  ! pressure, erg/cm3
         y0(i, idx_energy) = p1 / (gamma - 1d0) + y0(i, idx_rvx)**2 / y0(i, idx_rho) / 2.  ! update energy density, erg/cm3
      enddo

      ! advance time
      t = t + dt
      print *, t / tend
      if(t >= tend) exit
    end do

    ! copy final solutions to subroutine arguments
    rho = y0(:, idx_rho)
    vx = y0(:, idx_rvx) / rho
    p = (y0(:, idx_energy) - rho * vx**2 / 2.) * (gamma - 1d0)
    tgas = p / rho / kboltzmann * mu * pmass
    xall = y0(:, nhydro+1:nvars)

  end subroutine solve

end module hydro

! ****************
! driver code
program euler
  use hydro
  use krome_main
  use krome_user
  implicit none
  real*8,parameter::spy=3600.*24.*365  ! s / yr
  real*8::rho(ncell), vx(ncell), tgas(ncell), xall(ncell, krome_nmols), tend
  real*8::rhoL, rhoR, TgasL, TgasR
  integer::i

  ! density left and right conditions, g/cm3
  rhoL = 1d-20
  rhoR = 0.125 * rhoL
  ! temperature left and right conditions (to have pressure factor 0.1), K
  TgasL = 1.4d3
  TgasR = 0.1 * rhoL * TgasL / rhoR

  ! gas initially at rest, cm/s
  vx(:) = 0d0

  ! final integration time, s
  tend = 2.45d4 * spy

  ! set initial conditions
  rho(1:ncell/2) = rhoL
  rho(ncell/2+1:ncell) = RhoR
  tgas(1:ncell/2) = TgasL
  tgas(ncell/2+1:ncell) = TgasR

  ! initialize KROME
  call krome_init()

  ! set initial chemistry conditions
  xall = 1d-40
  xall(:, krome_idx_H) = 0.75615 * rho
  xall(:, krome_idx_E) = 4.4983d-8 * rho
  xall(:, krome_idx_Hj) = 8.1967d-5 * rho
  xall(:, krome_idx_He) = 0.24375 * rho
  xall(:, krome_idx_H2) = 1.5123d-6 * rho

  ! do the shock test
  call solve(rho, vx, tgas, xall, tend)

  ! write output to file
  do i=1,ncell
    write(22, '(99E17.8e3)') i * dx / pc2cm, rho(i), vx(i), tgas(i), &
      xall(i, krome_idx_H2), xall(i, krome_idx_Hk), xall(i, krome_idx_E)
  end do

  print *, "plot with gnuplot"
  print *, "load 'plot.gps'"
  print *, "done, bye!"

end program euler
