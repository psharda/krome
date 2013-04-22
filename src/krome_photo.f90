module krome_photo
contains
  subroutine krome_init_photo()
    use krome_commons

    print *,"Initializaing photoreactions..."
#KROME_photo_init_zero
#KROME_photo_qromos
#KROME_photo_heating_qromos

#KROME_photo_heating_print
    print *,"done!"

  end subroutine krome_init_photo

#KROME_photo_functions
#KROME_photo_heating_functions


  !********************
  function Jflux(nrg)
    !UV flux in eV/s/cm2/sr/Hz as a function of energy in eV
    implicit none
    real*8::Jflux,nrg
    Jflux = 6.2415d-10 * (13.6d0/nrg) / 4.d1 !6.241d-10eV = 1d-21erg
  end function Jflux

  !*******************************
  function intf(nrg)
    use krome_constants
    real*8::intf,nrg
    intf = 4.d0 * pi * Jflux(nrg) / nrg / planck_eV * eV_to_erg !* exp(-tau) 
  end function intf

  !********************
  function sigma_v96(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
    real*8::sigma_v96,energy_eV,sigma_0,Fy,yw,x,y,E0
    real*8::y0,y1,ya,P
    x = energy_eV/E0 - y0
    y = sqrt(x**2 + y1**2)
    Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
         * (1.d0+sqrt(y/ya))**(-P)
    sigma_v96 = 1d-18 * sigma_0 * Fy !cm2
  end function sigma_v96

  !********************
  function heat_v96(energy_eV,Eth,E0,sigma_0,ya,P,yw,y0,y1)
    use krome_constants
    real*8::heat_v96,energy_eV,sigma_0,Fy,yw,x,y,E0,Eth
    real*8::y0,y1,ya,P
    x = energy_eV/E0 - y0
    y = sqrt(x**2 + y1**2)
    Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
         * (1.d0+sqrt(y/ya))**(-P)
    heat_v96 = 1d-18 * sigma_0 * Fy * (energy_eV - Eth) !cm2*eV
  end function heat_v96

  !********************************
  subroutine qromos(func, sigma, a, b, ss, choose)
    !modified qromo to have sigma*function
    implicit none

    !     Parameters
    integer, parameter :: jmax = 18, jmaxp = jmax+1, km = 4, k = km+1
    real*8, parameter :: eps = 1.0e-06

    !     Arguments
    real*8 :: a, b, ss

    !     Externals
    real*8, external :: func,sigma
    external :: choose

    !     Locals
    real*8 :: s(jmaxp), h(jmaxp)
    real*8 :: dss
    integer :: j

    h(1) = 1.0

    do j = 1, jmax
       call choose(func,sigma,a,b,s(j),j)
       if (j.ge.k) then
          call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
          if (abs(dss) .lt. (eps*abs(ss)) ) return
       endif

       s(j+1) = s(j)
       h(j+1) = h(j)/9.0

       !     this is where the assumption of step tripling and an
       !     even error series is used.

    end do
    stop 'ERROR: too many steps in qromos.'

  end subroutine qromos

  !****************************************
  subroutine midsqls(funk,sigma,aa,bb,s,n)
    !     this routine is an exact replacement for midpnt except that it
    !     allows for an inverse square-root singularity int the integrand
    !     at the lower limit aa
    implicit none

    !     Arguments
    real*8 :: aa, bb, ss, s
    integer :: n

    real*8,external :: funk,sigma

    !     Locals
    real*8 :: a, b
    real*8 :: del, ddel, tnm, sum, x
    integer :: j
    integer :: it

    !     Statement function
    real*8 :: func

    func(x) = 2.0*x*funk(aa+x**2)*sigma(aa+x**2)

    save it

    b = sqrt(bb-aa)
    a = 0.0

    if (n.eq.1) then
       s = (b-a)*func(0.5*(a+b))
       it = 1
       !     2*it points will be added on the next refinement
    else
       tnm = it
       del = (b-a)/(3.0*tnm)
       ddel = del+del
       x = a+0.5*del
       sum = 0.0
       do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
       end do
       s = (s+(b-a)*sum/tnm)/3.0
       !     the new sum id combined with the old integral to give
       !     a refined integral.
       it = 3*it
    endif

    return
  end subroutine midsqls


  !********************************************

  subroutine qromo(func,a,b,ss,choose)

    !     Romberg integration on an open interval. returns as ss
    !     the integral of the function func from a to b. using
    !     any specified integrating sunroutine "choose" and
    !     Romberg's' method. normally "choose" will be an open
    !     formula, not evaluating the function at the endpoints.
    !     it is assumed that "choose" triples the number of steps
    !     on each call, and that its error series contains only
    !     even powers of the number of steps. the routines
    !     "midpnt", "midinf", "midsql", "midsqu", are possible
    !     choices for "choose".

    implicit none

    !     Parameters
    integer, parameter :: jmax = 18, jmaxp = jmax+1, km = 4, k = km+1
    real*8, parameter :: eps = 1.0e-06

    !     Arguments
    real*8 :: a, b, ss

    !     Externals
    real*8, external :: func
    external :: choose

    !     Locals
    real*8 :: s(jmaxp), h(jmaxp)
    real*8 :: dss
    integer :: j

    h(1) = 1.0

    do j = 1, jmax
       call choose(func,a,b,s(j),j)
       if (j.ge.k) then
          call polint(h(j-km),s(j-km),k,0.d0,ss,dss)
          if (abs(dss) .lt. (eps*abs(ss)) ) return
       endif
       s(j+1) = s(j)
       h(j+1) = h(j)/9.0
       !     this is where the assumption of step tripling and an
       !     even error series is used.
    end do
    stop 'too many steps.'
  end subroutine qromo

  !**********************************************
  subroutine polint(xa,ya,n,x,y,dy)

    !     this is a routine for polynomial interpolation or 
    !     extrapolation. geven arrays xa and ya, each of lenth
    !     n, and value x, this routine will return a value y,
    !     and an error estimator dy. if p(x) is the polynomial
    !     of degree n-1 such that p(xa_j) = ya_j, j = 1,...,n, then
    !     the returned value y = p(x).

    implicit none
    !     Parameter
    !     change nmax as desired to be the largest anticipated
    !     value of n.
    integer, parameter :: nmax = 10

    !     Arguments
    integer :: n
    real*8 :: xa(n), ya(n)
    real*8 :: x, y, dy

    !     Locals
    real*8 :: c(nmax), d(nmax)
    real*8 :: dif, dift, ho, hp, w, den
    integer :: ns
    integer :: i, m

    ns = 1
    dif = abs(x-xa(1))

    do i = 1, n
       !     here we find the index ns of the closest table entry.
       dift = abs(x-xa(i))
       if (dift.lt.dif) then
          ns = i
          dif = dift
       endif
       c(i) = ya(i)
       !     and initialize the tableau of c's' and d's'.
       d(i) = ya(i)
    end do
    y = ya(ns)
    !     this is the initial appromation to y.
    ns = ns-1
    do m = 1, n-1
       !     for each column of the tableau,
       do i = 1, n-m
          !     we loop over the current c's' and d's' and update them.
          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp
          if(den == 0.0) stop 'den  =  0'
          !     this error can only occur if two input xa's' are identical
          den = w/den
          d(i) = hp*den
          c(i) = ho*den
       end do

       if ((2*ns) < (n-m)) then
          !     after each column in the tableau is completed, we
          !     decided which correction , c or d, we want to add 
          !     to our accumulating value of y, i.e. which path to
          !     take through the tableau --- forking up or down.
          !     we do this in such a way as to take the most 
          !     "straight line" route through the tableau to its apex,
          !     updating ns accordingly to keep track of where we are.
          !     this route keeps the partial approximationa centered
          !     (insofar as possible) on the target x. the last dy
          !     added is thus the error indication.
          dy = c(ns+1)
       else
          dy = d(ns)
          ns = ns-1
       endif
       y = y+dy
    end do

    return
  end subroutine polint

  !****************************************
  subroutine midinf(funk,aa,bb,s,n)
    !     this subroutine is an exact replacement for "midpnt",
    !     i.e. returns as "s" the n-th stage of refinement of the
    !     integral of "func" from aa to bb, except that the function
    !     is evaluated at evenly spaced points in 1/x rather than 
    !     in x. this allows the upper limit bb to be as large and
    !     positive as the computer allows, or the lower limit aa
    !     to be as large and negative, but not both. aa and bb must
    !     have the same sign.
    implicit none
    !     Arguments
    integer :: n
    real*8 :: aa, bb, s
    real*8 :: funk
    external :: funk
    !     Locals
    real*8 :: a, b
    real*8 :: del, ddel, x, sum, tnm
    integer :: j
    integer :: it
    !     Statement function
    real*8 :: func

    func(x) = funk(1.0/x)/x**2
    save it
    b = 1.0/aa
    a = 1.0/bb

    if (n == 1) then
       s = (b-a)*func(0.5*(a+b))
       it = 1
    else
       tnm = it
       del = (b-a)/(3.0*tnm)
       ddel = del+del
       x = a+0.5*del
       sum = 0.0
       do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
       end do
       s = (s+(b-a)*sum/tnm)/3.0
       it = 3*it
    endif

    return
  end subroutine midinf


  !****************************************
  subroutine midpnt(func,a,b,s,n)
    !     this routine computes the n-th stage of refinement of an 
    !     extended midpoint rule. "func" is input as the name of the 
    !     function to be integrated between limits a and b, also input.
    !     when called with n = 1, the routine returns as "s" the crudest
    !     estimate of integral. subsequent calls with n = 2,3,... (in that
    !     sequential order) will improve the accuracy of "s" by adding
    !     (2/3)*3^{n-1} additional interior points. "s" should not be
    !     modified between sequential calls.

    implicit none
    !     Arguments
    real*8 :: a, b, s
    integer :: n
    real*8 :: func
    external :: func
    !     Locals
    real*8 :: del, ddel, tnm, sum, x
    integer :: j
    integer :: it

    save it

    if (n.eq.1) then
       s = (b-a)*func(0.5*(a+b))
       it = 1
       !     2*it points will be added on the next refinement
    else
       tnm = it
       del = (b-a)/(3.0*tnm)
       ddel = del+del
       x = a+0.5*del
       sum = 0.0
       do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
       end do
       s = (s+(b-a)*sum/tnm)/3.0
       !     the new sum id combined with the old integral to give
       !     a refined integral.
       it = 3*it
    endif

    return
  end subroutine midpnt

  !****************************************
  subroutine midsql(funk,aa,bb,s,n)
    !     this routine is an exact replacement for midpnt except that it
    !     allows for an inverse square-root singularity int the integrand
    !     at the lower limit aa
    implicit none
    !     Arguments
    real*8 :: aa, bb, ss, s
    integer :: n
    real*8 :: funk
    external :: funk
    !     Locals
    real*8 :: a, b
    real*8 :: del, ddel, tnm, sum, x
    integer :: j
    integer :: it
    !     Statement function
    real*8 :: func
    func(x) = 2.0*x*funk(aa+x**2)

    save it

    b = sqrt(bb-aa)
    a = 0.0
    if (n.eq.1) then
       s = (b-a)*func(0.5*(a+b))
       it = 1
       !     2*it points will be added on the next refinement
    else
       tnm = it
       del = (b-a)/(3.0*tnm)
       ddel = del+del
       x = a+0.5*del
       sum = 0.0
       do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
       end do
       s = (s+(b-a)*sum/tnm)/3.0
       !     the new sum id combined with the old integral to give
       !     a refined integral.
       it = 3*it
    endif

    return
  end subroutine midsql

  !****************************************
  subroutine midsqu(funk,aa,bb,s,n)
    !     this routine is an exact replacement for midpnt except that it
    !     allows for an inverse square-root singularity int the integrand
    !    at the upper limit bb
    implicit none
    !     Arguments
    real*8 :: aa, bb, s
    integer :: n

    real*8 :: funk
    external :: funk
    !     Locals
    real*8 :: a, b
    real*8 :: del, ddel, tnm, sum, x
    integer :: j
    integer :: it
    !     Statement function
    real*8 :: func

    func(x) = 2.0*x*funk(bb-x**2)

    save it

    b = sqrt(bb-aa)
    a = 0.0

    if (n.eq.1) then
       s = (b-a)*func(0.5*(a+b))
       it = 1
       !     2*it points will be added on the next refinement
    else
       tnm = it
       del = (b-a)/(3.0*tnm)
       ddel = del+del
       x = a+0.5*del
       sum = 0.0

       do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
       end do
       s = (s+(b-a)*sum/tnm)/3.0
       !     the new sum id combined with the old integral to give
       !     a refined integral.
       it = 3*it
    endif

    return
  end subroutine midsqu

end module krome_photo
