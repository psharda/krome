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
  real*8::flux(krome_nPhotoBins),Av

  call krome_init()

  Tgas(:) = 7d2 !gas temperature in K (all cells)
  spy = krome_seconds_per_year !s
  pc = 3.085d18 !1 pc in cm
  Rsun = 7d10 !sun radius in cm
  rmin = pc*1d-3 !min radius: cm
  rmax = pc*1d1 !box radius: cm
  Rstar = 66.1*Rsun !emitting sphere radius
  ntot = 1d4

  !set cosmic rays
  call krome_set_user_crate(1.3d-17)

  !init abundances
  xall(:,:) = 0d0
  xall(:,krome_idx_H) = ntot
  !xall(:,krome_idx_H2) = ntot/3.*2.
  xall(:,krome_idx_C) = 2.5d-4*ntot
  xall(:,krome_idx_O) = 4.7d-4*ntot
  xall(:,krome_idx_He) = 8.5d-2*ntot

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
  call krome_set_photoBin_draineLog(6d0, 12d0)
  call krome_load_opacity_table("opdust.dat",unitEnergy="micron")
  !scale the draine flux by 10 times
  call krome_photoBin_scale(1.75d4)
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
     print *,t/spy
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
        Av = Av + dr * sum(xall(i,:))/1.8d21
        !opacity from photochemitry
        op(:) = op(:)*exp(-krome_get_opacity_size_d2g(x(:),Tgas(i), dr,1d-2))
        !geometrical flux attenuation
        opx(:) = op(:)/4d0/krome_pi!/4d0*(rstar/r)**2
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
        av = av + (r-rold)*sum(x) / 1.8d21
        write(66,'(99E17.8e3)') t/spy,av,Tgas(i),x(:)
        rold = r
     end do
     write(66,*)

     !exit after tend
     if(t>1d9*spy) exit
  end do

  !loop in grid points to write final conditions
  Av = 0d0
  rold = 0d0
  do i=1,ngrid
     x(:) = xall(i,:)
     r = rx(i)
     av = av + (r-rold)*sum(x) / 1.8d21
     write(55,'(99E17.8e3)') t/spy,Av,Tgas(i),x(:)/ntot
     rold = r
  end do

  print *,"done!"

end program main
