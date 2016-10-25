!#################################################################
!This test consists of a gas with an initial population of dust,
! both silicon and carbon based, with 10 size-bins each.
!The gas-to-dust ratio is 10e-5 and the dust distribution
! follows an MRN power-law. Tgas evolves (see main loop).
!#################################################################
program test_krome

  use krome_main
  use krome_user

  integer,parameter::nd=krome_ndust,imax=100
  real*8::x(krome_nmols),Tgas,t,dt,spy,xH,tend,vgas,xi(krome_nmols)
  real*8::xdust(nd),adust(nd),xdusti(nd),data(imax,nd),dataT(imax)
  integer::i,j

  spy = 365.*24.*3600. !seconds per year
  Tgas = 1d1 !gas temperature (K)
  xH = 1d6 !Hydrogen density

  !initialize krome
  call krome_init()

  x(:) = 1d-20 !default species abundance (cm-3)

  !number densities (cm-3)
  x(KROME_idx_H)     = 1d0 * xH  !H
  x(KROME_idx_Hj)    = 1d-4 * xH   !H+
  x(KROME_idx_H2)    = 1d-5 *xH   !H2
  x(KROME_idx_C)     = 1d-4 *xH    !C
  x(KROME_idx_Si)    = 1d-4 *xH    !Si

  !compute electrons (globally neutral)
  x(KROME_idx_E) = krome_get_electrons(x(:))

  xi(:) = x(:) !store the initial amount of species
  dt = 1d3*spy !time-step (s)

  !output file header
  write(66,*) "#Tgas mC mSi"

  !loop on Tgas
  do i=1,imax
     x(:) = xi(:)
     !re-init dust
     call krome_init_dust_distribution(x(:),1d-2)
     Tgas = 1d1**((i-1)*(4d0-1d0) / (imax-1) + 1d0)
     call krome_set_Tdust(Tgas) !Tdust is coupled with Tgas
     if(mod(i,10)==0) print *,Tgas !print every 10

     !call KROME
     call krome(x(:),Tgas,dt)

     !get evolved dust size
     adust(:) = krome_get_dust_size()
     !compute mass density in the dust bins
     data(i,:) = 4d0*krome_pi/3d0*adust(:)**3*2.2d0*krome_get_dust_distribution()
     !write mass density in the gas phase, C and Si
     write(66,*) Tgas, x(krome_idx_C)*12d0*krome_p_mass, &
          x(krome_idx_Si)*28d0*krome_p_mass
     write(88,*) Tgas,sum(data(i,:)) + x(krome_idx_C)*12d0*krome_p_mass + &
          x(krome_idx_Si)*28d0*krome_p_mass
     !store Tgas
     dataT(i) = Tgas
   end do

   !write dust mass density evolution
   write(77,*) "#dustType Tgas rhodust"
   !loop on dust bins
   do j=1,nd
      !loop on data
      do i=1,imax
         !dust type, Tgas, dust mass density
         write(77,*) (j-1)/(nd/2),dataT(i),data(i,j)
      end do
      write(77,*)
   end do

   print *,"load 'plot.gps' in gnuplot to show the results"
   print *,"That's all! have a nice day!"

end program test_krome
