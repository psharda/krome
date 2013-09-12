module mods
contains

  subroutine lin(xmin,xmax,arr)
    real*8::arr(:),xmin,xmax
    integer::i
    do i=1,size(arr)
       arr(i) = i*(xmax-xmin)/size(arr) + xmin
    end do
  end subroutine lin


  !***************
  subroutine rhs ( x_num, x, t, value )

    !*****************************************************************80
    !
    !! RHS_TEST01 evaluates the right hand side for problem 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    25 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) X_NUM, the number of nodes.
    !
    !    Input, real ( kind = 8 ) X(X_NUM), the node coordinates.
    !
    !    Input, real ( kind = 8 ) T, the current time.
    !
    !    Output, real ( kind = 8 ) VALUE(X_NUM), the source term.
    !
    implicit none

    integer ( kind = 4 ) x_num

    real ( kind = 8 ) t
    real ( kind = 8 ) value(x_num)
    real ( kind = 8 ) x(x_num)

    value(1:x_num) = 0.0D+00

    return
  end subroutine rhs

  !************************
  subroutine fbc ( x_num, x, t, h )

    !*****************************************************************80
    !
    !! BC_TEST01 evaluates the boundary conditions for problem 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    25 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) X_NUM, the number of nodes.
    !
    !    Input, real ( kind = 8 ) X(X_NUM), the node coordinates.
    !
    !    Input, real ( kind = 8 ) T, the current time.
    !
    !    Input, real ( kind = 8 ) H(X_NUM), the current heat values.
    !
    !    Output, real ( kind = 8 ) H(X_NUM), the current heat values, after boundary
    !    conditions have been imposed.
    !
    implicit none

    integer ( kind = 4 ) x_num

    real ( kind = 8 ) h(x_num)
    real ( kind = 8 ) t
    real ( kind = 8 ) x(x_num)

    h(1)  = 0.d0
    h(x_num) = 0.d0

    return
  end subroutine fbc

end module mods
