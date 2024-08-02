module krome_fit
contains

  !*****************************
  subroutine init_anytab3D(filename,x,y,z,f,xmul,ymul,zmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8::rout(4)
    integer::i,j,k,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(f,1)) then
       print *,"ERROR: in init_anytab3D x size differs from f(x,y,z)"
       stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(f,2)) then
       print *,"ERROR: in init_anytab3D y size differs from f(x,y,z)"
       stop
    end if

    !check the size of the Z input array
    if(size(z).ne.size(f,3)) then
       print *,"ERROR: in init_anytab3D z size differs from f(x,y,z)"
       stop
    end if

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: in init_anytab3D file ",trim(filename)," not found!"
       stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
       read(unit,'(a)') row_string
       if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
       print *,"ERROR: file "//filename//" should"
       print *," contain the number of grid points"
       print *," per dimension in the format"
       print *,"  XX, YY, ZZ"
       print *,row_string
       stop
    end if

    !loop to read file (3rd dimension of f() is
    ! first in the tables. i.e. tables are z,x,y,
    ! while f() is x,y,z
    do i=1,size(z)
       do j=1,size(x)
          do k=1,size(y)
             read(unit,*,iostat=ios) rout(:)
             y(k) = rout(3)
             f(j,k,i) = rout(4)
          end do
          x(j) = rout(2)
          read(unit,*,iostat=ios) !skip blanks
       end do
       z(i) = rout(1)
       read(unit,*,iostat=ios) !skip blanks
       if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))
    zmul = 1d0/(z(2)-z(1))

  end subroutine init_anytab3D

  !********************************************
  !load 2d tables from filename
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:,:),xmul,ymul
    real*8::rout(3)
    integer::i,j,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
       print *,"ERROR: in init_anytab2D x size differs from z"
       stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
       print *,"ERROR: in init_anytab2D y size differs from z"
       stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
       stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
       read(unit,'(a)') row_string
       if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
       print *,"ERROR: file "//filename//" should"
       print *," contain the number of rows and "
       print *," columns in the format"
       print *,"  RR, CC"
       print *,row_string
       stop
    end if

    !loop to read file
    do i=1,size(x)
       do j=1,size(y)
          read(unit,*,iostat=ios) rout(:)
          y(j) = rout(2)
          z(i,j) = rout(3)
       end do
       x(i) = rout(1)
       read(unit,*,iostat=ios) !skip blanks
       if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !********************************************
  !load 1d tables from filename
  subroutine init_anytab1D(filename,x,y,xmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),xmul
    real*8::rout(2)
    integer::i,ios,unit

    !check the size of the X input array
    if(size(x) /= size(y)) then
       print *,"ERROR: in init_anytab1D x size differs from y"
       stop
    end if

    if (krome_mpi_rank <= 1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios /= 0) then
       print *,"ERROR: in init_anytab1D file ",trim(filename)," not found!"
       stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
       read(unit,'(a)') row_string
       if(row_string(1:1)/="#") exit
    end do

    ! !check if first line is OK
    ! if(scan(row_string,",")==0) then
    !    print *,"ERROR: file "//filename//" should"
    !    print *," contain the number of rows and "
    !    print *," columns in the format"
    !    print *,"  RR, CC"
    !    print *,row_string
    !    stop
    ! end if

    !loop to read file
    do i=1,size(x)
      read(unit,*,iostat=ios) rout(:)
      y(i) = rout(2)
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios /= 0) exit

    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))

  end subroutine init_anytab1D

  !******************************
  !test 2d fit and save to file
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
       do j=1,size(y)
          xx = x(i)
          yy = y(j)
          zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
          write(unit1,*) xx,yy,zz
          write(unit2,*) x(i),y(j),z(i,j)
       end do
       write(unit1,*)
       write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_anytab2D

  !******************************
  !test 2d interpolation and save to file
  subroutine test_interpolate2D(fname,x,y,z)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:)
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
       do j=1,size(y)
          xx = x(i)
          yy = y(j)
          zz = interpolate2D(x(:),y(:),z(:,:),xx,yy)
          write(unit1,*) xx,yy,zz
          write(unit2,*) x(i),y(j),z(i,j)
       end do
       write(unit1,*)
       write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_interpolate2D

  !*****************************
  function fit_anytab3D(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D

  !******************************
  !return 2d fit at xx,yy
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    implicit none
    real*8::fit_anytab2D
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  !return 1d fit at xx
  function fit_anytab1D(x,z,xmul,xx)
    real*8,intent(in)::x(:),z(:),xmul,xx
    real*8::fit_anytab1D,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  function fit_anytab3D_linlinlog(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D_linlinlog,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D_linlog(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D_linlog(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D_linlinlog = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D_linlinlog

  !***************************
  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
       z1 = log10(-z1)
       z2 = log10(-z2)
       fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
       return
    end if

    if(z1>0d0 .and. z2>0d0) then
       z1 = log10(z1)
       z2 = log10(z2)
       fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
       return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !1D interpolation at x0 for x(:) in z(:)
  !Added by Piyush Sharda in 2024 for CIE cooling from Gnat and Ferland 2012
  function interpolate1D(x, z, x0)
    real*8 :: x(:), z(:)
    real*8 :: x0
    real*8 :: interpolate1D
    integer :: i
    real*8 :: t

    ! Find index i such that x(i) <= x0 < x(i+1)
    i = 1
    do while (i < size(x) - 1 .and. x0 > x(i + 1))
      i = i + 1
    end do

    ! Compute interpolation weight
    t = (x0 - x(i)) / (x(i + 1) - x(i))

    ! Perform linear interpolation
    interpolate1D = (1 - t) * z(i) + t * z(i + 1)
    
  end function interpolate1D

  !*****************************
  !2D interpolation at (x0,y0) for (x(:), y(:)) in z(:,:)
  !Added by Piyush Sharda in 2024 for CO shielding
  function interpolate2D(x, y, z, x0, y0)
    real*8 :: interpolate2D
    real*8 :: x(:), y(:), z(size(x), size(y))
    real*8 :: x0, y0
    real*8 :: f
    integer :: i, j
    real*8 :: t, u

    ! Find indices i and j such that x(i) <= x0 < x(i+1) and y(j) <= y0 < y(j+1)
    i = 1
    do while (i < size(x) - 1 .and. x0 > x(i + 1))
      i = i + 1
    end do

    j = 1
    do while (j < size(y) - 1 .and. y0 > y(j + 1))
      j = j + 1
    end do

    ! Compute interpolation weights
    t = (x0 - x(i)) / (x(i + 1) - x(i))
    u = (y0 - y(j)) / (y(j + 1) - y(j))

    ! Perform bilinear interpolation
    interpolate2D = (1 - t) * (1 - u) * z(i, j) + t * (1 - u) * z(i + 1, j) + &
        (1 - t) * u * z(i, j + 1) + t * u * z(i + 1, j + 1)
        
  end function interpolate2D

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
       b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
       c(1) = 0d0
       d(1) = 0d0
       b(2) = b(1)
       c(2) = 0d0
       d(2) = 0d0
       return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
       d(i) = x(i+1) - x(i)
       b(i) = 2d0*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
       ispline = y(1)
       return
    end if

    if(u>=x(n)) then
       ispline = y(n)
       return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
       k = (i+j)/2
       if(u<x(k)) then
          j=k
       else
          i=k
       end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_fit
