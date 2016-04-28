!In this test we have a semi-infinite slab
! in a box with size rmin-rmax defined by ngrid grid points.
!The light is attenuated via opacity
! determined by the photochemical xsecs present in the network
! and dust opacity loaded from file
!The emission is 10 times the draine flux in the energy
! interval [5,13.6] eV
program main
  use krome_main
  use krome_user
  use krome_user_commons
  implicit none
  integer,parameter::ngrid=300 !grid points
  integer::i,j
  real*8::Tgas(ngrid),dt,x(krome_nmols),spy,pc,Rsun,ntot
  real*8::xall(ngrid,krome_nmols),rmin,rmax,r,t,rx(ngrid)
  real*8::Rstar,Tstar,dr,rold,op(krome_nPhotoBins)
  real*8::opx(krome_nPhotoBins),energy(krome_nPhotoBins)
  real*8::flux(krome_nPhotoBins),Av,fav,tend

  !fav = 1.8d21
  !Av/NH conversion factor
  fav = 1d0/6.289d-22

  call krome_init()

  Tgas(:) = 5d1 !gas temperature in K (all cells)
  ntot = 1d3 !total density
  spy = krome_seconds_per_year !s
  pc = 3.085d18 !1 pc in cm
  Rsun = 7d10 !sun radius in cm
  rmin = fav/ntot*1d-4 !min radius: cm
  rmax = fav/ntot*3d0 !box radius: cm
  Rstar = 66.1*Rsun !emitting sphere radius
  tend = 1d6*spy

  !set cosmic rays
  call krome_set_user_crate(5d-17)

  !init abundances
  xall(:,:) = 0d0
  xall(:,krome_idx_H) = ntot/3.
  xall(:,krome_idx_H2) = ntot/3.*2.
  xall(:,krome_idx_C) = 1d-4*ntot !2.5d-4*ntot
  xall(:,krome_idx_O) = 3d-4*ntot !4.7d-4*ntot
  xall(:,krome_idx_He) = 1d-1*ntot !8.5d-2*ntot

  !initial time step
  dt = .1*spy
  t = 0d0

  !loop on the grid points
  do i=1,ngrid
     !compute the position of the grid points
     rx(i) = 1d1**((i-1)*(log10(rmax)-log10(rmin))/(ngrid-1) &
          + log10(rmin))
  end do

  !uses a draine flux
  call krome_set_photoBin_draineLin(5d0, 13.6d0)
  !load energy (micron), opacity (cm2/g) table from file
  call krome_load_opacity_table("opacityDraineR35.dat",unitEnergy="micron")
  !scale the draine flux by 10 times
  call krome_photoBin_scale(1d1)
  call krome_photoBin_store()
  !store energy value for plot
  energy(:) = krome_get_photoBinE_mid()

  !loop on time
  do
     !open a file to write the opacity (replaced every time step)
     open(45,file="opacity.dat",status="replace")
     !increase time-step size
     dt = dt * 1.2
     !increase total time
     t = t + dt
     !print time
     print *,t/tend*1d2,"%"
     rold = 0d0
     !initial opacity is 1 (thin)
     op(:) = 1d0
     Av = 0d0
     !loop on grid points
     do i=1,ngrid
        !restore the flux to original every frist grid point
        call krome_photoBin_restore()
        r = rx(i)
        !grid spacing
        dr = r-rold
        !compute Av
        Av = Av + dr * sum(xall(i,:))/fav
        !opacity from photochemitry
        op(:) = op(:)*exp(-krome_get_opacity_size_d2g(x(:),Tgas(i), dr,1d-2))
        !semi-infinite slab (krome assumes 4*pi)
        opx(:) = op(:) / 4d0
        !store flux for plot
        flux(:) = krome_get_photoBinJ()
        !loop on photobins to write energy-dependent flux
        do j=1,krome_nPhotoBins
           write(45,*) Av,energy(j),op(j)*flux(j)
        end do
        write(45,*)
        !scale the flux according to opacity
        call krome_photoBin_scale_array(opx(:))
        !call KROME to do chemistry
        x(:) = xall(i,:)
        call krome(x(:),Tgas(i),dt)
        xall(i,:) = x(:)
        rold = r
     end do
     close(45)

     !write the time evolution for the ngrid/2 cell
     x(:) = xall(ngrid/2,:)
     write(77,'(99E17.8e3)') t/spy,Tgas(ngrid/2),x(:)

     rold = 0d0
     av = 0d0
     !write full evolution (all cells, all time-steps)
     do i=1,ngrid
        x(:) = xall(i,:)
        r = rx(i)
        av = av + (r-rold)*sum(x) / fav
        write(66,'(99E17.8e3)') t/spy,av,Tgas(i),x(:)
        rold = r
     end do
     write(66,*)

     !exit after tend
     if(t>tend) exit
  end do

  !loop in grid points to write final conditions
  Av = 0d0
  rold = 0d0
  do i=1,ngrid
     x(:) = xall(i,:)
     r = rx(i)
     av = av + (r-rold)*sum(x) / fav
     write(55,'(99E17.8e3)') t/spy,Av,Tgas(i),x(:)
     rold = r
  end do

  print *,"done!"

end program main
