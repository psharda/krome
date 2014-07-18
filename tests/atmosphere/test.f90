!######################################################################
!This test couple the standard diffusion equation for a gas
! with the solution of the ODE system for chemical species
! based on the network of Kasting & Donahue 1980.
!The initial conditions are taken from Segura et al. 2003.
!In order to keep the chemical time-scale similar to the 
! dynamical time a normalization of the chemical species
! abundances is applied. For further details refer to the KROME paper.
!######################################################################
program test
  use krome_main
  use krome_user
  implicit none
  integer,parameter::nmax=64,nread=59
  real*8::r(nmax), tt, dt, x(krome_nmols),h(nmax)
  real*8::rout(2),rout2(4),tgas(nmax),n(nmax,krome_nmols),n1(nmax)
  real*8::datar(nread),ntot,dtin,tmax,xx(krome_nmols),xmax
  integer::i,j,iout,istep,jmax


  call krome_init()

  tmax = 2d4 !total simulation time

  !read eddy values from file
  open(33,file="eddy.dat",status="old")
  !loop on layers
  do i=1,nmax
     read(33,*) rout(:)
     r(i) = rout(2) !second column
  end do
  close(33)

  r(:) = r(:) / maxval(r) * 5d-1 !normalize to 0.5

  !read initial layers data
  open(33,file="layers_data",status="old")
  read(33,*) !read header
  !loop on layers
  do i=1,nmax
     read(33,*) iout, rout2(:)
     h(i) = rout2(1) * 1d5 !km->cm
     tgas(i) = rout2(3) !temperature in K
  end do
  close(33)

  !read species
  open(33,file="init_spec.dat",status="old")
  do i=1,nmax
     x(:) = 0.d0 !default
     read(33,*) datar(:) !read file line
     x(krome_idx_H2O) = datar(1)
     x(krome_idx_O_1D) = datar(2)
     x(krome_idx_OH) = datar(3)
     x(krome_idx_H) = datar(4)
     x(krome_idx_H2) = datar(5)
     x(krome_idx_O) = datar(6)
     x(krome_idx_O2) = datar(7)
     x(krome_idx_O3) = datar(8)
     x(krome_idx_HO2) = datar(9)
     x(krome_idx_H2O2) = datar(11)
     x(krome_idx_N2) = datar(12)
     x(krome_idx_O_3P) = datar(13)
     x(krome_idx_CO) = datar(15)
     x(krome_idx_CO2) = datar(16)
     x(krome_idx_HCO) = datar(17)
     x(krome_idx_H2CO) = datar(18)
     x(krome_idx_CH4) = datar(19)
     x(krome_idx_CH3OOH) = datar(20)
     x(krome_idx_H3CO) = datar(21)
     x(krome_idx_N2O) = datar(22)
     x(krome_idx_NO) = datar(24)
     x(krome_idx_HNO3) = datar(25)
     x(krome_idx_NO2) = datar(26)
     x(krome_idx_N) = datar(27)
     x(krome_idx_CH3) = datar(28)
     x(krome_idx_CH2_1) = datar(30)
     x(krome_idx_CH2_3) = datar(31)
     x(krome_idx_CH3O2) = datar(32)
     x(krome_idx_NO3) = datar(33)
     !copy to layer
     n(i,:) = x(:)
  end do
  close(33)

  n(:,:) = n(:,:) / maxval(n) !normalize
  print *,"Wait, it takes a while..."
  
  dt = 1d0 !minval(r, MASK = r>0) / 1d3
  tt = 0 !absolute time
  istep = 0 !count steps
  do
     !do diffusion
     do j=1,krome_nmols
        ntot = sum(n(:,j)) 
        do i = 2, nmax - 1
           n1(i) = n(i,j) + r(i) * (n(i - 1,j) - 2 * n(i,j) + n(i + 1,j))
        end do
        n1(1) = n1(2)
        n1(nmax) = n1(nmax-1)
        n(:,j) = n1(:) 
     end do

     !do chemistry
     !$omp parallel do private(i, dtin, x) schedule(dynamic,2)
     do i=1,nmax
        x(:) = n(i,:)
        dtin = dt
        call krome(x(:), Tgas(i), dtin)
        n(i,:) = x(:)
     end do
     
     !dump every 1000 steps
     if(mod(istep,1000)==0) then
        print '(F11.2,a2)',tt/tmax*1d2," %"
        do i = 1,nmax
           x(:) = n(i,:)
           write(55,'(999E17.8e3)') tt,h(i),x(:)
        end do
        write(55,*)
     end if
     tt = tt + dt !increase time
     if(tt>tmax) exit !exit condition
     istep = istep + 1 !increase timestep
  end do

  print *,"done!"
  print *,"type:"
  print *," >load 'plot.gps'"
  print *,"in gnuplot to show the results"

end program test
