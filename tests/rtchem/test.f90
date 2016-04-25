!In this test we have an emmiting sphere that cast rays
! in a box with size rmin-rmax defined by ngrid grid points.
!The light is attenuated via geometric decay and opacity
! determined by the photochemical xsecs present in the network
!The emission is 10 times the draine flux in the energy
! interval [11.6,13.6] eV
program main
  use krome_main
  use krome_user
  implicit none
  integer,parameter::ngrid=30 !grid points
  integer::i,j
  real*8::Tgas(ngrid),dt,x(krome_nmols),spy,pc,Rsun,ntot
  real*8::xall(ngrid,krome_nmols),rmin,rmax,r,t,rx(ngrid)
  real*8::Rstar,Tstar,dr,rold,op(krome_nPhotoBins)
  real*8::opx(krome_nPhotoBins),energy(krome_nPhotoBins)
  real*8::flux(krome_nPhotoBins)

  call krome_init()

  Tgas(:) = 5d1 !gas temperature in K (all cells)
  spy = krome_seconds_per_year !s
  pc = 3.085d18 !1 pc in cm
  Rsun = 7d10 !sun radius in cm
  rmin = pc*1d-5 !min radius: cm
  rmax = pc*0.005 !box radius: cm
  Rstar = 66.1*Rsun !emitting sphere radius
  ntot = 1d3

  !set cosmic rays
  call krome_set_user_crate(1.3d-17)

  !init abundances
  xall(:,:) = 0d0
  xall(:,krome_idx_H) = ntot/3.
  xall(:,krome_idx_H2) = ntot/3.*2.
  xall(:,krome_idx_C) = 1d-4
  xall(:,krome_idx_O) = 3d-4

  !initial time step
  dt = .1*spy
  t = 0d0

  !loop on the grid points
  do i=1,ngrid
     !compute the position of the grid points
     rx(i) = 1d1**((i-1)*(log10(rmax)-log10(rmin))/(ngrid-1) &
          + log10(rmin))
  end do

  !loop on time
  do
     !open a file to write the opacity (replaced every time step)
     open(45,file="opacity.dat",status="replace")
     !uses a draine flux
     call krome_set_photoBin_draineLog(11.6d0, 13.6d0)
     !scale the draine flux by 10 times
     call krome_photoBin_scale(1d1)
     !store energy value for plot
     energy(:) = krome_get_photoBinE_mid()
     !increase time-step size
     dt = dt * 1.2
     !increase total time
     t = t + dt
     print *,t/spy
     rold = 0d0
     !initial opacity is 1 (thin)
     op(:) = 1d0
     !loop on grid points
     do i=1,ngrid
        !restore the flux to original every frist grid point
        call krome_photoBin_restore()
        r = rx(i)
        !grid spacing
        dr = r-rold
        !opacity from photochemitry
        op(:) = op(:)*exp(-krome_get_opacity_size(x(:),Tgas(i), dr))
        !geometrical flux attenuation
        opx(:) = op(:)/4d0*(rstar/r)**2
        !store flux for plot
        flux(:) = krome_get_photoBinJ()
        !loop on photobins to write energy-dependent flux
        do j=1,krome_nPhotoBins
           write(45,*) r/pc,energy(j),op(j)*flux(j)
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

     !write full evolution (all cells, all time-steps)
     do i=1,ngrid
        x(:) = xall(i,:)
        r = rx(i)
        write(66,'(99E17.8e3)') t/spy,r/pc,Tgas(i),x(:)
     end do
     write(66,*)

     !exit after tend
     if(t>1d9*spy) exit
  end do

  !loop in grid points to write final conditions
  do i=1,ngrid
     x(:) = xall(i,:)
     r = rx(i)
     write(55,'(99E17.8e3)') t/spy,r/pc,Tgas(i),x(:)
  end do

  print *,"done!"

end program main
