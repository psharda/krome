!###################################################

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use krome_constants
  implicit none
  real*8::Tgas,x(krome_nmols),dt,spy,ntot,t
  integer::i,imax,j,jmax

  spy = krome_seconds_per_year
  call krome_init() !init krome (mandatory)

  ntot = 1d5

  imax = 10
  jmax = 10
  print '(a8,a12)',"%","ntot"
  write(44,*) "#ntot Tgas cooling"
  do j=1,jmax
     ntot = 1d1**((j-1)*(9d0-3d0)/(jmax-1)+3d0)
     print '(F8.1,E12.3)',j*1d2/jmax,ntot
     do i=1,imax
        Tgas = 1d1**((i-1)*(4d0-1d0)/(imax-1)+1d0)

        x(:) = 1.d-40

        x(krome_idx_H)  = ntot          !H
        x(krome_idx_H2) = 1d-6*ntot    !H2
        x(krome_idx_Hj) = 1d-4*ntot    !H+
        x(krome_idx_HE) = 0.0775d0*ntot !He

        !set cosmic rays
        call krome_set_user_crate(1.3d-17)

        !rescale metallicity for neutral metals (C,Fe,Si,O)
        call krome_scale_Z(x(:), 1d-1)

        call krome_set_user_Av((ntot*1d-3)**(2./3.))

        x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
        x(krome_idx_C)  = 1d-40

        x(krome_idx_e) = krome_get_electrons(x(:))

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
        write(44,'(99E17.8e3)') ntot,Tgas,krome_coolingO2(x(:),Tgas)/ntot
        write(55,'(99E17.8e3)') Tgas, x(:)
     end do
     write(44,*)
  end do

  print *,"Test OK!"

end program test

