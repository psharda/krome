program test
  use krome_main
  use krome_user
  implicit none
  integer,parameter::nmax=64,nread=59
  real*8::r(nmax), tt, dt, x(krome_nmols)
  real*8::xmax,xmin,rout(2),rout2(4),tgas(nmax),n(nmax,krome_nmols),n1(nmax)
  real*8::datar(nread),ntot,dtin
  integer::i,j,iout,istep


  call krome_init()

  !read eddy values
  open(33,file="eddy.dat",status="old")
  do i=1,nmax
     read(33,*) rout(:)
     r(i) = rout(2) 
  end do
  close(33)

  r(:) = r(:) / maxval(r) * 1d-1

  !read initial layers data
  open(33,file="layers_data",status="old")
  read(33,*)
  do i=1,nmax
     read(33,*) iout, rout2(:)
     !x(i) = rout2(1) * 1d5
     tgas(i) = rout2(3)
  end do
  close(33)

  !read species
  open(33,file="init_spec.dat",status="old")
  do i=1,nmax
     x(:) = 0.d0
     read(33,*) datar(:)
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

     x(:) = 0d0
     
     x(krome_idx_OH) = i * 1d0 / nmax
     n(i,:) = x(:)
  end do
  close(33)

  n(:,:) = n(:,:) / maxval(n)
  print *,maxval(n),minval(n)
  
  xmax = maxval(x)
  xmin = minval(x)

  dt = 1d-1 !minval(r, MASK = r>0) / 1d3
  tt = 0
  istep = 0
  do
     do j=1,krome_nmols
        ntot = sum(n(:,j))
        do i = 2, nmax - 1
           n1(i) = n(i,j) + r(i) * (n(i - 1,j) - 2 * n(i,j) + n(i + 1,j))
        end do
        n1(1) = n1(2)
        n1(nmax) = n1(nmax-1)
        n(:,j) = n1(:) !/ sum(n1(:)) * ntot
     end do

     do i=1,nmax
        x(:) = n(i,:)
        dtin = dt
        !call krome(x(:), Tgas(i), dtin)
        n(i,:) = x(:)
     end do

     if(mod(istep,1000)==0) then
        do i = 1,nmax
           write(55,*) tt,i,n(i,krome_idx_H2O)
        end do
        write(55,*)
     end if
     tt = tt + dt
     if(istep>int(1e6)) exit
     istep = istep + 1
  end do

end program test
