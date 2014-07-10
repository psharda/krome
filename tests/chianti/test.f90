!###################################################
!This is a test to plot cooling function using e.g.
! the information from the Chianti database.
! http://www.chiantidatabase.org
! It computes a cooling curve similar to the one 
! shown in Gnat+Ferland 2011 ApJS 199 20 (fig.2)

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  real*8::Tgas,dt,x(krome_nmols),rho,spy,t,j21s(3),j21
  integer::j,k,kmax

  spy = krome_seconds_per_year !use shorter variable for this constant

  j21s(:) = (/0d0,1d0,1d2/) !list of j21 values
  kmax = 30  !number of temperature intervals

  call krome_init() !init krome (mandatory)

  !loop on j21 values
  do j=1,size(j21s)

     !init radiation field between 1d1 and 1d2 using 1/E profile
     call krome_set_photoBin_J21log(1d1, 1d2)

     j21 = j21s(j) !select J21 value
     call krome_photoBin_scale(j21) !scale radiation according to j21
     print *,"running with j21=",j21

     !loop on temeprature values
     do k=1,kmax

        !init abundances
        x(:) = 0d0  !default abundances
        x(krome_idx_Cj) = 1d0
        x(krome_idx_e) = 1d0!krome_get_electrons(x(:))

        Tgas = 1d1**((k-1)*(8.-4.)/(kmax-1)+4.) !gas temperature (K)
        if(mod(k,10)==0) print '(F7.2,a2)',(k-1)*1d2/(kmax-1),"%"

        dt = spy !initial time-step (s)
        t = 0d0 !start time (s)

        !switch off thermo to get equilibirum at constant temperature
        call krome_thermo_off()
        
        !equilibrium
        do
           dt = dt * 1.1d0 !increase time-step
           t = t + dt !advance time
           !electron conservation
           x(krome_idx_e) = krome_get_electrons(x(:))
           call krome(x(:), Tgas, dt) !call KROME
           write(66,'(I5,99E17.8)') j, t/spy, Tgas, x(:) !dump
           if(t>1d8*spy) exit !exit when 1d8 years reached
        end do

        !write the abundances of the species and corresponding cooling
        write(77,'(I5,999E17.8e3)') j,Tgas,x(:)
        write(55,'(I5,999E17.8e3)') j,Tgas,krome_get_cooling(x,Tgas)
        write(66,*)
     end do
     write(55,*)
     write(77,*)
     write(66,*)
  end do

  print *,"done!"
  print *," for a plot in gnuplot type"
  print *," load 'plot.gps'"

end program test

