module KROME_coolingGH
#IFKROME_useCoolingGH
    ! ------------------------------------------------------------
    !
    ! This module: Cooling and Heating Functions, table
    ! reader for Gnedin and Hollon cooling tables
    ! Language:  Fortran 77
    !
    !  UPDATED by Troels Haugboelle to Fortran 90, explicit kinds,
    !  and encapsulated in a module
    !
    !  Copyright (c) 2012 Nick Gnedin
    !  All rights reserved.
    !
    !  Redistribution and use in source and binary forms, with or without
    !  modification, are permitted provided that the following conditions
    !  are met:
    !
    !  Redistributions of source code must retain the above copyright
    !  notice, this list of conditions and the following disclaimer.
    !
    !  Redistributions in binary form must reproduce the above copyright
    !  notice, this list of conditions and the following disclaimer in the
    !  documentation and/or other materials provided with the distribution.
    !
    !  Neither the name of Nick Gnedin nor the names of any contributors
    !  may be used to endorse or promote products derived from this software
    !  without specific prior written permission.
    !
    !  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    !  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    !  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    !  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
    !  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    !  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    !  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    !  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
    !  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    !  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    !  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    ! ------------------------------------------------------------
    private
    !  Table dimensions
    integer,parameter::NT=81, NX=13, NP1=24, NP2=21, NP3=16, ND=3789
    !  Number of components per T-D bin
    integer,parameter:: NC=6, NICH=12, NRCH=13
    !  Mode of table lookup
    integer::mode
    !  Boundaries
    integer::np(3)
    !  Data and index blocks
    integer:: indx(NP1,NP2,NP3)
    real*4::data(NC,NT,NX,ND)
    ! Indices
    real*8::altval(NT), altmin, altstp, &
         xval(NX), xmin, xmax, xstp, &
         qmin(3), qmax(3), qstp(3)
    public :: frtInitCF, frtCFCache, frtCFGetLn, frtGetCF
    public :: NICH, NRCH
  contains

    subroutine frtInitCF(m,fname)
      implicit none
      integer                      :: m
      character(len=*), intent(in) :: fname
      !  Internally used unit number
      integer, parameter :: IOCF=97
      integer :: i, j, k, id, ix, it, ic, lt, ld, lp1, lp2, lp3, lp4, lx
      real*4 :: q1, q2
      real*4 :: altval4(NT), xmin4, xmax4, qmin4(3), qmax4(3)

      mode = m

      open(unit=IOCF, file=fname, status='old', form='unformatted', err=100)
      read(IOCF,err=100) lt, ld, lp1, lp2, lp3, lp4, &
           (qmin4(j),j=1,3), q1, (qmax4(j),j=1,3), q2, lx, xmin4, xmax4

      if(lt.ne.NT .or. ld.ne.ND .or. lx.ne.NX .or. lp1.ne.NP1 .or. &
           lp2.ne.NP2 .or. lp3.ne.NP3 .or. lp4.ne.1 .or. ld.eq.0) then
         write(0,*) 'RT::InitCF: fatal error, corrupted table:'
         write(0,*) '> NT= in file: ', lt, ' in code: ', NT
         write(0,*) '> NX= in file: ', lx, ' in code: ', NX
         write(0,*) '> ND= in file: ', ld, ' in code: ', ND
         write(0,*) '> NP1= in file: ', lp1, ' in code: ', NP1
         write(0,*) '> NP2= in file: ', lp2, ' in code: ', NP2
         write(0,*) '> NP3= in file: ', lp3, ' in code: ', NP3
         write(0,*) '> NP4= in file: ', lp4, ' in code: ', 1
         close(IOCF)
         m = -1
         stop
      end if

      qmin=qmin4; qmax=qmax4; xmin=xmin4; xmax=xmax4

      np(1) = lp1
      np(2) = lp2
      np(3) = lp3

      do i=1,3
         if(np(i) .gt. 1) then
            qstp(i) = (qmax(i)-qmin(i))/(np(i)-1)
         else
            qstp(i) = 1.0
         end if
      end do

      xstp = (xmax-xmin)/(NX-1)
      do i=1,NX
         xval(i) = xmin + xstp*(i-1)
      end do

      read(IOCF,err=100) (altval4(i),i=1,NT)
      altval = altval4
      !  Internally use natural log
      do i=1,NT
         altval(i) = altval(i)*log(10.0)
      end do
      altmin = altval(1)
      altstp = altval(2) - altval(1)

      read(IOCF,err=100) (((indx(i,j,k),i=1,lp1),j=1,lp2),k=1,lp3)

      do id=1,ld
         read(IOCF,err=100) (((data(ic,it,ix,id),ic=1,NC),it=1,NT),ix=1,NX)
      end do

      do id=1,ND
         do ix=1,NX
            do it=1,NT
               do ic=1,NC
                  data(ic,it,ix,id) = log(1d-37+abs(data(ic,it,ix,id)))
               end do
            end do
         end do
      end do

      close(IOCF)

      if(.false.) then
         write(6,*) 'RT::InitCF: Table size = ', &
              NC*(NT*NX*ND/256/1024) + (NP1*NP2*NP3/256/1024), ' MB'
      endif

      m = 0

      return

100   m = -1

    end subroutine frtInitCF
    !
    !  Decode the interpolated function
    !
    !#define IXL      ich(3)
    !#define IXU      ich(4)
    !#define IPP(j)   ich(4+j)
    !#define WXL      rch(6)
    !#define WXU      rch(7)
    !#define WPL(j)   rch(7+j)
    !#define WPS(j)   rch(10+j)
    !
    subroutine frtCFPick(it,ich,rch,cfun,hfun)
      implicit none
      integer :: it
      integer :: ich(:)
      real*8 :: rch(:)
      real*8 :: cfun, hfun
      !
      integer :: ic, j
      real*8 :: v(NC), q(8), a0, a1, a2, Z
      !
      do ic=1,NC
         do j=1,8
            q(j) = rch(6)*data(ic,it,ich(3),ich(4+j)) + &
                 rch(7)*data(ic,it,ich(4),ich(4+j))
         end do
         v(ic) = exp( &
              rch(10)*(rch(9)*(rch(8)*q( 1)+rch(11)*q( 2))+   &
              rch(12)*(rch(8)*q( 3)+rch(11)*q( 4))) + &
              rch(13)*(rch(9)*(rch(8)*q( 5)+rch(11)*q( 6))+   &
              rch(12)*(rch(8)*q( 7)+rch(11)*q( 8))))

      end do

      a0 = v(1)
      a1 = v(2)
      a2 = v(3)
      v(2) = 2*a1 - 0.5*a2 - 1.5*a0
      v(3) = 0.5*(a0+a2) - a1

      a0 = v(4)
      a1 = v(5)
      a2 = v(6)
      v(5) = 2*a1 - 0.5*a2 - 1.5*a0
      v(6) = 0.5*(a0+a2) - a1

      Z = rch(5)

      if(mode .eq. 1) then
         cfun = (Z*v(3)+v(2))*Z
         hfun = (Z*v(6)+v(5))*Z
      else
         cfun = (Z*v(3)+v(2))*Z + v(1)
         hfun = (Z*v(6)+v(5))*Z + v(4)
      end if

    end subroutine frtCFPick

    !***********************
    !  Cache some table information into arrays iCache and rCache
    subroutine frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
      implicit none
      real*8, intent(in)   :: den, Z, Plw, Ph1, Pg1, Pc6
      real*8, dimension(:) :: rch
      integer, dimension(:) :: ich
      integer :: ierr
      !
      real*8  :: q(3), qh1, qg1, qc6, dl, w
      integer :: j, il(3), is(3)

      !  Convert from nb to nH from Cloudy models
      dl = max(1.0e-10,den*(1.0-0.02*Z)/1.4)
      ierr = 0

      if(Plw .gt. 0.0) then
         qh1 = log10(1.0e-37+Ph1/Plw)
         qg1 = log10(1.0e-37+Pg1/Plw)
         qc6 = log10(1.0e-37+Pc6/Plw)

         q(1) = log10(1.0e-37+Plw/dl)
         q(2) = 0.263*qc6 + 0.353*qh1 + 0.923*qg1
         q(3) = 0.976*qc6 - 0.103*qh1 - 0.375*qg1

         !  qmin, qstp, etc are boundaries of cells, not their centers
         do j=1,3
            w = 0.5 + (q(j)-qmin(j))/qstp(j)
            il(j) = int(w) + 1
            if(w .gt. il(j)-0.5) then
               is(j) = il(j) + 1
            else
               is(j) = il(j) - 1
            endif
            rch(10+j) = abs(il(j)-0.5-w)
            rch(7+j) = 1 - rch(10+j)

            if(np(j) .gt. 1) then
               if(max(il(j),is(j)) .gt. np(j)) ierr =  j
               if(min(il(j),is(j)) .lt.     1) ierr = -j
            endif

            if(il(j) .lt. 1) il(j) = 1
            if(is(j) .lt. 1) is(j) = 1
            if(il(j) .gt. np(j)) il(j) = np(j)
            if(is(j) .gt. np(j)) is(j) = np(j)
         enddo

      else

         ierr = -1
         do j=1,3
            il(j) = 1
            is(j) = 1
            rch(7+j) = 1
            rch(10+j) = 0
         enddo

      endif

      !  Density interpolation is still CIC
      w = (log10(dl)-xmin)/xstp
      ich(3) = int(w) + 1
      if(ich(3) .lt.  1) ich(3) = 1
      if(ich(3) .ge. NX) ich(3) = NX-1
      ich(4) = ich(3) + 1
      rch(6) = max(0.0,min(1.0,ich(3)-w))
      rch(7) = 1.0 - rch(6)

      !  Do not forget C-to-F77 index conversion
      ich(5) = 1 + indx(il(1),il(2),il(3))
      ich(6) = 1 + indx(is(1),il(2),il(3))
      ich(7) = 1 + indx(il(1),is(2),il(3))
      ich(8) = 1 + indx(is(1),is(2),il(3))
      ich(9) = 1 + indx(il(1),il(2),is(3))
      ich(10) = 1 + indx(is(1),il(2),is(3))
      ich(11) = 1 + indx(il(1),is(2),is(3))
      ich(12) = 1 + indx(is(1),is(2),is(3))

      rch(5) = Z

      !  Clear temperature cache
      ich(1) = 0
      ich(2) = 0
    end subroutine frtCFCache

    ! Get the cooling and heating functions for given
    !  ln(T) from the cached dat
    subroutine frtCFGetLn(alt,ich,rch,cfun,hfun)
      real*8 :: alt
      real*8,  dimension(:) :: rch
      integer, dimension(:) :: ich
      real*8 :: cfun, hfun
      real*8, dimension(NC) :: v
      integer      :: il, iu
      real*8 :: ql, qu
      !
      il = int((alt-altmin)/altstp*0.99999) + 1
      if(il .lt.  1) il = 1
      if(il .ge. NT) il = NT-1
      iu = il + 1
      ql = max(0.0,min(1.0,(altval(iu)-alt)/altstp))
      qu = 1.0 - ql
      !
      !  Shift cache lines as needed
      !
      if(ich(1) .eq. iu) then
         ich(1) = 0
         ich(2) = iu
         rch(3) = rch(1)
         rch(4) = rch(2)
      endif
      if(ich(2) .eq. il) then
         ich(1) = il
         ich(2) = 0
         rch(1) = rch(3)
         rch(2) = rch(4)
      endif

      !  Update the cache
      if(ich(1) .ne. il) then
         ich(1) = il
         call frtCFPick(il,ich,rch,cfun,hfun)
         rch(1) = cfun
         rch(2) = hfun
      endif
      if(ich(2) .ne. iu) then
         ich(2) = iu
         call frtCFPick(iu,ich,rch,cfun,hfun)
         rch(3) = cfun
         rch(4) = hfun
      endif

      cfun = ql*rch(1) + qu*rch(3)
      hfun = ql*rch(2) + qu*rch(4)
    end subroutine frtCFGetLn

    !  Get the cooling and heating functions for T
    subroutine frtGetCF(tem,den,Z,Plw,Ph1,Pg1,Pc6,cfun,hfun,ierr)
      implicit none
      real*8 :: tem, den, Z, Plw, Ph1, Pg1, Pc6, cfun, hfun
      integer :: ierr
      ! Cache arrays
      real*8, dimension(NRCH) :: rch
      integer,      dimension(NICH) :: ich

      call frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
      call frtCFGetLn(log(max(1.0,tem)),ich,rch,cfun,hfun)

      if (ierr > 0) then
        print *,'frtGetCF: Problems with caching Gnedin and Hollon cooling table.'
        print *,'frtGetCF:  Stopping. Error code :', ierr
        print *,'frtGetCF: Input variables ntot, Tgas, PLW, PHI, PHeI, PCVI: ', den, tem, Plw, Ph1, Pg1, Pc6
        stop
      end if

    end subroutine frtGetCF

#ENDIFKROME
  end module KROME_coolingGH

