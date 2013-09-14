program test
  use libheat
  use mods
  !external rhs,bc
  integer,parameter::np=10,tnum=200
  real*8::x(np),t(tnum),dt,cc,h(np),hout(np)
  real*8::k,tmin,tmax,xmin,xmax,cfl, htot
  integer::i,j
  tmin = 0.d0
  tmax = 5.d1
  k = 0.02d0
  xmin = 0.d0
  xmax = 1.d0
  call lin(xmin,xmax,x(:))
  call lin(tmin,tmax,t(:))

  
   
  dt = (tmax-tmin)/tnum
  h(1) = .0d0
  h(np) = .0d0
  do j=2,np-1
     h(j) = -4.d0*x(j)**2 + 4.d0*x(j)
  end do

  htot = sum(h)

  call fcfl(k, tnum, tmin, tmax, np, xmin, xmax, cc )
  do i=1,tnum
     call solver(np,x(:),t(i),dt,cc,rhs,fbc,h(:),hout(:))
     hout(np) = hout(np-1)
     hout(1) = hout(2) 
     h(:) = hout(:) / sum(hout) * htot
     print *,sum(h)
     do j=1,np
        write(22,'(99E12.3)') x(j),t(i),h(j) 
     end do
     write(22,*)
  end do

  print *,"bznznzbbzgggg"
end program test
