!##################################################
! Test based on [R2008] Reinhardt, V. et al.
! JPhysChemA 112, 1712 (2008)
! Slow manifold trajectories example
! This program evolves the chemical network for
! jmax*kmax initial conditions of H2O and H2
!##################################################
program test_krome
  use krome_main
  use krome_subs
  use krome_user_commons
  use krome_user

  integer,parameter::nspec=6
  integer::i,imax,j,k,jmax,kmax
  real*8::x(nspec),Tgas,rho,dt,xold(nspec),t,xi(nspec)
  character*50::names(nspec)

  !grid points
  jmax = 14
  kmax = 7

  Tgas = 1.2d3 !gas Temperature (K), dummy
  print *,"solving..."

  write(66,*) "#time "//krome_get_names_header()
  write(77,*) "#time "//krome_get_names_header()

  !loop on H2O points
  do j=1,jmax
     !loop on H2 points
     do k=1,kmax
        !set initial conditions
        xi(krome_idx_OH)  = 0d0 !OH
        xi(krome_idx_O)   = 0d0 !O
        xi(krome_idx_H2O) = (j-1) * (0.65 - 0.05) / (jmax-1) + 0.05 !H2O
        xi(krome_idx_H2)  = (k-1) * (0.9-0.3) / (kmax-1) + 0.3 !H2
        !O2 and H from Eqn.(16) of [R2008]
        xi(krome_idx_O2)  = (1d0 - xi(krome_idx_H2o) - xi(krome_idx_OH) &
             - xi(krome_idx_O)) / 2d0
        xi(krome_idx_H)   = 2d0 -xi(krome_idx_OH) -2d0*xi(krome_idx_H2O) &
             -2d0*xi(krome_idx_H2)
        !use only initial values on the lower part of the H2-H2O space
        if(xi(krome_idx_H2)>.96d0-xi(krome_idx_H2O)) cycle
        dt = 1d-10 !initial time-step (s)
        t = 0d0 !initial time (s)
        x(:) = xi(:) !init abundances
        write(77,'(999E12.3e3)') t,x(:) !dump initial conditions to file
        !loop on time
        do
           xold(:) = x(:) !store fractions
           dt = dt * 1.2d0 !increase time-step
           call krome(x(:), Tgas, dt) !call KROME
           t = t + dt !increase time
           write(66,'(999E12.3e3)') t,x(:) !dump evolution
           if(sum((xold(:)-x(:))**2)/nspec<1d-20) exit !check the steady-state
        end do
        write(66,*) !blank line
        write(66,*) !blank line
     end do
  end do

  !print final abundances to screen
  print *,"Final equilibrium abundances"
  names(:) = krome_get_names()
  do j=1,nspec
     print '(a4,E11.3)',names(j),x(j)
  end do

  print *,"done!"
  print *,"Output written in fort.66"
  print *,"Default gnuplot: load 'plot.gps'"

  print *,"that's all! have a nice day!"

end program test_krome
