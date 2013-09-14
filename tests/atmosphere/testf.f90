program test

  integer, parameter :: NMAX = 64,nspec=3
  real*8::T(NMAX), T1(NMAX), r(nmax), tt, dt, x(nmax)
  real*8::xmax,xmin,rout(2),rout2(4),tgas(nmax),n(nmax,nspec),n1(nmax)
  real*8::nout(nspec)
  integer::i,j,iout,istep


  open(33,file="eddy.dat",status="old")
  do i=1,nmax
     read(33,*) rout(:)
     r(i) = rout(2) * 1d-7
  end do
  close(33)

  open(33,file="layers_data",status="old")
  read(33,*)
  do i=1,nmax
     read(33,*) iout, rout2(:)
     x(i) = rout2(1) * 1d5
     tgas(i) = rout2(3)
  end do
  close(33)

  open(33,file="init_spec.dat",status="old")
  do i=1,nmax
     read(33,*) nout(:)
     !print *,nout(:)
     n(i,:) = nout(:)
  end do
  close(33)

  xmax = maxval(x)
  xmin = minval(x)


  dt = minval(r) / 1d6
  tt = 0
  istep = 0
  do
     do j=1,nspec
        do i = 2, nmax - 1
           !print *,n(i,j)
           n1(i) = n(i,j) + r(i) * (n(i - 1,j) - 2 * n(i,j) + n(i + 1,j))
        end do
        n1(1) = n1(2)
        n1(nmax) = n1(nmax-1)
        n(:,j) = n1(:)
     end do


     if(mod(istep,100)==0) then
        print *,sum(n1)
        do i = 1,nmax
           write(55,*) tt,i,n(i,1)
        end do
        write(55,*)
     end if
     tt = tt + dt
     istep = istep + 1
  end do

end program test
