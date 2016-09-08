!###################################################
!This is a test to plot cooling function using e.g.
! the information from the Chianti database.
! http://www.chiantidatabase.org
! It computes a cooling curve similar to the one
! shown in Gnat+Ferland 2012 ApJS 199 20 (Fig.2)
! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1111.6980

program test
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  implicit none
  real*8::Tgas,dt,x(krome_nmols),spy,t,j21s(3),j21
  integer, parameter :: kmax=30, imaxx=5000
  integer::i,j,k,imax(kmax)
  real*8::res66(2+krome_nmols,imaxx,kmax), res55(2,kmax)
  real*8::res77(1+krome_nmols,kmax) !temp arrs for output

  spy = krome_seconds_per_year !use shorter variable for this constant

  j21s(:) = (/0d0,1d0,1d2/) !list of j21 values

  call krome_init() !init krome (mandatory)

  !loop on j21 values
  do j=1,size(j21s)

     !$omp parallel private(j21)

     !init radiation field between 1d1 and 1d2 using 1/E profile
     call krome_set_photoBin_J21log(1d1, 1d2)

     j21 = j21s(j) !select J21 value
     call krome_photoBin_scale(j21) !scale radiation according to j21

     !$omp end parallel

     print *,"running with j21=",j21

     !loop on temperature values
     !$omp parallel do private(i,k,x,Tgas,dt,t) schedule(dynamic,1)
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
        i=0
        do
           i=i+1
           dt = dt * 1.1d0 !increase time-step
           t = t + dt !advance time
           !electron conservation
           x(krome_idx_e) = krome_get_electrons(x(:))
           call krome(x(:), Tgas, dt) !call KROME
           if (i <= imaxx) then
             res66(:,i,k) = (/ t/spy, Tgas, x(:) /)
             imax(k) = i
           endif
           if(t>1d8*spy) exit !exit when 1d8 years reached
        end do

        !write the abundances of the species and corresponding cooling
        res77(:,k) = (/ Tgas, x(:) /)
        res55(:,k) = (/ Tgas, krome_get_cooling(x,Tgas) /)
     end do

     !dump results on a file
     write(55,'(a)') "#J21 Tgas cooling"
     do k=1,kmax
        do i=1, imax(k)
           write(66,'(I5,99E17.8)') j, res66(:,i,k)
        enddo
        write(77,'(I5,999E17.8e3)') j, res77(:,k)
        write(55,'(999E17.8e3)') j21s(j), res55(:,k)
        write(66,*)
     enddo
     write(55,*)
     write(77,*)
     write(66,*)
  end do

  print *,"done!"
  print *," for a plot in gnuplot type"
  print *," load 'plot.gps'"

end program test

