  module KROME_cooling
#KROME_header
    integer,parameter::coolTab_n=int(1e2)
#KROME_nZrate
    real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
    real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
#KROME_escape_vars
#KROME_coolingZ_popvars
  contains

    !*******************
    function cooling(n,inTgas)
      use krome_commons
      implicit none
      real*8::n(:),inTgas,cooling,Tgas

      Tgas = inTgas
      cooling = sum(get_cooling_array(n(:),Tgas))

    end function cooling

    !*******************************
    function get_cooling_array(n, Tgas)
      use krome_commons
      implicit none
      real*8::n(:), Tgas
      real*8::get_cooling_array(ncools),cools(ncools)
      real*8::f1,f2,smooth

      f1 = 1d0
      f2 = 1d0

      !returns cooling in erg/cm3/s
      cools(:) = 0d0

#IFKROME_useCoolingH2
      cools(idx_cool_H2) = cooling_H2(n(:), Tgas) #KROME_floorH2
#ENDIFKROME

#IFKROME_useCoolingH2GP
      cools(idx_cool_H2GP) = cooling_H2GP(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingHD
      cools(idx_cool_HD) = cooling_HD(n(:), Tgas) #KROME_floorHD
#ENDIFKROME

#IFKROME_useCoolingdH
      cools(idx_cool_dH) = cooling_dH(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDust
      cools(idx_cool_dust) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDustNoTdust
      cools(idx_cool_dust) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDustTabs
      cools(idx_cool_dust) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingExpansion
      cools(idx_cool_exp) = cooling_expansion(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCO
      cools(idx_cool_CO) = cooling_CO(n(:), Tgas) #KROME_floorCO
#ENDIFKROME

#IFKROME_useCoolingOH
      cools(idx_cool_OH) = cooling_OH(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingH2O
      cools(idx_cool_H2O) = cooling_H2O(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingHCN
      cools(idx_cool_HCN) = cooling_HCN(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingAtomic
      cools(idx_cool_atomic) = cooling_Atomic(n(:), Tgas) #KROME_floorAtomic
#ENDIFKROME

#IFKROME_useCoolingCompton
      cools(idx_cool_compton) = cooling_compton(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingFF
      cools(idx_cool_ff) = cooling_ff(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingContinuum
      cools(idx_cool_cont) = cooling_continuum(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingGH
      !this parameter controls the smoothness of the
      ! merge between the two cooling functions
      smooth = 1./500.

      !smoothing functions | f1+f2=1
      f1 = (tanh(smooth*(Tgas-2.5d3))+1.d0)*0.5d0
      f2 = (tanh(smooth*(-Tgas+2.5d3))+1.d0)*0.5d0

      ! Only compute cooling terms if it really matters
      if (f2 > 1d-20) then
#ENDIFKROME

#IFKROME_useCoolingZ
      cools(idx_cool_Z) = f2 * ( cooling_Z(n(:), Tgas) #KROME_floorZ )
#ENDIFKROME

#IFKROME_useCoolingCIE
      cools(idx_cool_CIE) = f2 * cooling_CIE(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingZCIE
      cools(idx_cool_ZCIE) = f2 * ( cooling_Z_CIE(n(:), Tgas) #KROME_floorZ_CIE )
#ENDIFKROME

#IFKROME_useCoolingZCIENOUV
      cools(idx_cool_ZCIENOUV) = f2 * ( cooling_Z_CIENOUV(n(:), Tgas) #KROME_floorZ_CIENOUV )
#ENDIFKROME

#IFKROME_useCoolingZExtended
      !floor is inside the function
      cools(idx_cool_ZExtend) = f2 * cooling_ZExtended(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingGH
      endif

      if (f1 > 1d-20) cools(idx_cool_GH) = f1 * cooling_GH(n(:), Tgas)
#ENDIFKROME

      cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

      get_cooling_array(:) = cools(:)

    end function get_cooling_array

#IFKROME_useCoolingCO
    !***************************
    ! CO cooling: courtesy of K.Omukai (Nov2014)
    ! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
    ! see also Omukai+2010 (bit.ly/1HIaGcn)
    ! H and H2 collisions
    function cooling_CO(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      implicit none
      integer,parameter::imax=coolCOn1
      integer,parameter::jmax=coolCOn2
      integer,parameter::kmax=coolCOn3
      integer::i,j,k
      real*8,parameter::eps=1d-5
      real*8::cooling_CO,n(:),inTgas
      real*8::v1,v2,v3,prev1,prev2,cH
      real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
      real*8::x1(imax),x2(jmax),x3(kmax)
      real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
      real*8::v1min,v1max,v2min,v2max,v3min,v3max

      !local copy of limits
      v1min = coolCOx1min
      v1max = coolCOx1max
      v2min = coolCOx2min
      v2max = coolCOx2max
      v3min = coolCOx3min
      v3max = coolCOx3max

      !local copy of variables arrays
      x1(:) = coolCOx1(:)
      x2(:) = coolCOx2(:)
      x3(:) = coolCOx3(:)

      ixd1(:) = coolCOixd1(:)
      ixd2(:) = coolCOixd2(:)
      ixd3(:) = coolCOixd3(:)

      !local variables
      v3 = num2col(n(idx_CO),n(:)) !CO column density
      cH = n(idx_H) + n(idx_H2)
      if (n(idx_H) < 0. .or. n(idx_H2) < 0.) then
        cH = n_global(idx_H) + n_global(idx_H2)
        ! Set red flag if not due to small excursion in n_H2
        ! Only relevant for low temperatures, bc otherwise H is ionised
        ! And CO cooling only relevant at high densities
        if (abs(n(idx_H2)) / (abs(n(idx_H)) + 1d-40) > 1d-6 .and. &
                                               inTgas < 5d4 .and. &
            abs(n(idx_H)) > abs(n(idx_Hj)) .and. sum(n(1:nmols)) > 100.) &
                 red_flag = ibset(red_flag,5)
      endif

      v2 = cH
      v1 = inTgas !Tgas

      !logs of variables
      v1 = log10(v1)
      v2 = log10(v2)
      v3 = log10(v3)

      !default value erg/s/cm3
      cooling_CO = 0d0

      !check limits
      if(v1>=v1max) v1 = v1max*(1d0-eps)
      if(v2>=v2max) v2 = v2max*(1d0-eps)
      if(v3>=v3max) v3 = v3max*(1d0-eps)

      if(v1<v1min) return
      if(v2<v2min) return
      if(v3<v3min) return

      !gets position of variable in the array
      i = (v1-v1min)*coolCOdvn1+1
      j = (v2-v2min)*coolCOdvn2+1
      k = (v3-v3min)*coolCOdvn3+1

      !precompute shared variables
      prev1 = (v1-x1(i))*ixd1(i)
      prev2 = (v2-x2(j))*ixd2(j)

      !linear interpolation on x1 for x2,x3
      vv1 = prev1 * (coolCOy(k,j,i+1) - &
           coolCOy(k,j,i)) + coolCOy(k,j,i)
      !linear interpolation on x1 for x2+dx2,x3
      vv2 = prev1 * (coolCOy(k,j+1,i+1) - &
           coolCOy(k,j+1,i)) + coolCOy(k,j+1,i)
      !linear interpolation on x2 for x3
      vv12 = prev2 * (vv2 - vv1) + vv1

      !linear interpolation on x1 for x2,x3+dx3
      vv3 = prev1 * (coolCOy(k+1,j,i+1) - &
           coolCOy(k+1,j,i)) + coolCOy(k+1,j,i)
      !linear interpolation on x1 for x2+dx2,x3+dx3
      vv4 = prev1 * (coolCOy(k+1,j+1,i+1) - &
           coolCOy(k+1,j+1,i)) + coolCOy(k+1,j+1,i)
      !linear interpolation on x2 for x3+dx3
      vv34 = prev2 * (vv4 - vv3) + vv3

      !linear interpolation on x3
      xLd = (v3-x3(k))*ixd3(k)*(vv34 - &
           vv12) + vv12

      !CO cooling in erg/s/cm3
      cooling_CO = 1d1**xLd * cH * n(idx_CO)

    end function cooling_CO

    !************************
    subroutine init_coolingCO()
      use krome_commons
      implicit none
      integer::ios,iout(3),i,unit
      real*8::rout(4)

      if (krome_mpi_rank<=1) print *,"load CO cooling..."
      open(newunit=unit,file="coolCO.dat",status="old",iostat=ios)
      !check if file exists
      if(ios.ne.0) then
         print *,"ERROR: problems loading coolCO.dat!"
         stop
      end if

      do
         read(unit,*,iostat=ios) iout(:),rout(:) !read line
         if(ios<0) exit !eof
         if(ios/=0) cycle !skip blanks
         coolCOx1(iout(1)) = rout(1)
         coolCOx2(iout(2)) = rout(2)
         coolCOx3(iout(3)) = rout(3)
         coolCOy(iout(3),iout(2),iout(1)) = rout(4)
      end do

      !store inverse of the differences
      ! to speed up interpolation
      do i=1,coolCOn1-1
         coolCOixd1(i) = 1d0 / (coolCOx1(i+1)-coolCOx1(i))
      end do
      do i=1,coolCOn2-1
         coolCOixd2(i) = 1d0 / (coolCOx2(i+1)-coolCOx2(i))
      end do
      do i=1,coolCOn3-1
         coolCOixd3(i) = 1d0 / (coolCOx3(i+1)-coolCOx3(i))
      end do

      coolCOx1min = minval(coolCOx1)
      coolCOx1max = maxval(coolCOx1)
      coolCOx2min = minval(coolCOx2)
      coolCOx2max = maxval(coolCOx2)
      coolCOx3min = minval(coolCOx3)
      coolCOx3max = maxval(coolCOx3)

      coolCOdvn1 = (coolCOn1-1) / (coolCOx1max-coolCOx1min)
      coolCOdvn2 = (coolCOn2-1) / (coolCOx2max-coolCOx2min)
      coolCOdvn3 = (coolCOn3-1) / (coolCOx3max-coolCOx3min)

    end subroutine init_coolingCO
#ENDIFKROME

#IFKROME_useCoolingOH
    !***************************
    ! OH cooling: courtesy of K.Omukai (Mar2019)
    ! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
    ! see also Omukai+2010 (bit.ly/1HIaGcn)
    ! H and H2 collisions
    function cooling_OH(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      implicit none
      integer,parameter::imax=coolOHn1
      integer,parameter::jmax=coolOHn2
      integer,parameter::kmax=coolOHn3
      integer::i,j,k
      real*8,parameter::eps=1d-5
      real*8::cooling_OH,n(:),inTgas
      real*8::v1,v2,v3,prev1,prev2,cH
      real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
      real*8::x1(imax),x2(jmax),x3(kmax)
      real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
      real*8::v1min,v1max,v2min,v2max,v3min,v3max

      !local copy of limits
      v1min = coolOHx1min
      v1max = coolOHx1max
      v2min = coolOHx2min
      v2max = coolOHx2max
      v3min = coolOHx3min
      v3max = coolOHx3max

      !local copy of variables arrays
      x1(:) = coolOHx1(:)
      x2(:) = coolOHx2(:)
      x3(:) = coolOHx3(:)

      ixd1(:) = coolOHixd1(:)
      ixd2(:) = coolOHixd2(:)
      ixd3(:) = coolOHixd3(:)

      !local variables
      v3 = num2col(n(idx_OH),n(:)) !OH column density
      cH = n(idx_H) + n(idx_H2)

      v2 = cH
      v1 = inTgas !Tgas

      !logs of variables
      v1 = log10(v1)
      v2 = log10(v2)
      v3 = log10(v3)

      !default value erg/s/cm3
      cooling_OH = 0d0

      !check limits
      if(v1 >= v1max) v1 = v1max * (1d0 - eps)
      if(v2 >= v2max) v2 = v2max * (1d0 - eps)
      if(v3 >= v3max) v3 = v3max * (1d0 - eps)

      if(v1 < v1min) return
      if(v2 < v2min) return
      if(v3 < v3min) return

      !gets position of variable in the array
      i = (v1 - v1min) * coolOHdvn1 + 1
      j = (v2 - v2min) * coolOHdvn2 + 1
      k = (v3 - v3min) * coolOHdvn3 + 1

      !precompute shared variables
      prev1 = (v1 - x1(i)) * ixd1(i)
      prev2 = (v2 - x2(j)) * ixd2(j)

      !linear interpolation on x1 for x2,x3
      vv1 = prev1 * (coolOHy(k, j, i+1) - &
           coolOHy(k,j,i)) + coolOHy(k,j,i)
      !linear interpolation on x1 for x2+dx2,x3
      vv2 = prev1 * (coolOHy(k, j+1, i+1) - &
           coolOHy(k, j+1, i)) + coolOHy(k, j+1, i)
      !linear interpolation on x2 for x3
      vv12 = prev2 * (vv2 - vv1) + vv1

      !linear interpolation on x1 for x2,x3+dx3
      vv3 = prev1 * (coolOHy(k+1, j, i+1) - &
           coolOHy(k+1, j, i)) + coolOHy(k+1, j, i)
      !linear interpolation on x1 for x2+dx2,x3+dx3
      vv4 = prev1 * (coolOHy(k+1, j+1, i+1) - &
           coolOHy(k+1,j+1,i)) + coolOHy(k+1, j+1, i)
      !linear interpolation on x2 for x3+dx3
      vv34 = prev2 * (vv4 - vv3) + vv3

      !linear interpolation on x3
      xLd = (v3 - x3(k)) * ixd3(k) * (vv34 - &
           vv12) + vv12

      !OH cooling in erg/s/cm3
      cooling_OH = 1d1**xLd * cH * n(idx_OH)

    end function cooling_OH

    ! ************************
    subroutine init_coolingOH()
      use krome_commons
      implicit none
      integer::ios,iout(3),i,unit
      real*8::rout(4)

      if (krome_mpi_rank<=1) print *,"load OH cooling..."
      open(newunit=unit, file="coolOH.dat", status="old", iostat=ios)
      ! check if file exists
      if(ios.ne.0) then
         print *, "ERROR: problems loading coolOH.dat!"
         stop
      end if

      do
         read(unit,*,iostat=ios) iout(:),rout(:) !read line
         if(ios<0) exit !eof
         if(ios/=0) cycle !skip blanks
         coolOHx1(iout(1)) = rout(1)
         coolOHx2(iout(2)) = rout(2)
         coolOHx3(iout(3)) = rout(3)
         coolOHy(iout(3), iout(2), iout(1)) = rout(4)
      end do

      !store inverse of the differences
      ! to speed up interpolation
      do i=1,coolOHn1-1
         coolOHixd1(i) = 1d0 / (coolOHx1(i+1) - coolOHx1(i))
      end do
      do i=1,coolOHn2-1
         coolOHixd2(i) = 1d0 / (coolOHx2(i+1) - coolOHx2(i))
      end do
      do i=1,coolOHn3-1
         coolOHixd3(i) = 1d0 / (coolOHx3(i+1) - coolOHx3(i))
      end do

      coolOHx1min = minval(coolOHx1)
      coolOHx1max = maxval(coolOHx1)
      coolOHx2min = minval(coolOHx2)
      coolOHx2max = maxval(coolOHx2)
      coolOHx3min = minval(coolOHx3)
      coolOHx3max = maxval(coolOHx3)

      coolOHdvn1 = (coolOHn1 - 1) / (coolOHx1max - coolOHx1min)
      coolOHdvn2 = (coolOHn2 - 1) / (coolOHx2max - coolOHx2min)
      coolOHdvn3 = (coolOHn3 - 1) / (coolOHx3max - coolOHx3min)

    end subroutine init_coolingOH
#ENDIFKROME

#IFKROME_useCoolingH2O
    !***************************
    ! H2O cooling: courtesy of K.Omukai (Mar2019)
    ! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
    ! see also Omukai+2010 (bit.ly/1HIaGcn)
    ! H and H2 collisions
    function cooling_H2O(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      implicit none
      integer,parameter::imax=coolH2On1
      integer,parameter::jmax=coolH2On2
      integer,parameter::kmax=coolH2On3
      integer::i,j,k
      real*8,parameter::eps=1d-5
      real*8::cooling_H2O,n(:),inTgas
      real*8::v1,v2,v3,prev1,prev2,cH
      real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
      real*8::x1(imax),x2(jmax),x3(kmax)
      real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
      real*8::v1min,v1max,v2min,v2max,v3min,v3max

      !local copy of limits
      v1min = coolH2Ox1min
      v1max = coolH2Ox1max
      v2min = coolH2Ox2min
      v2max = coolH2Ox2max
      v3min = coolH2Ox3min
      v3max = coolH2Ox3max

      !local copy of variables arrays
      x1(:) = coolH2Ox1(:)
      x2(:) = coolH2Ox2(:)
      x3(:) = coolH2Ox3(:)

      ixd1(:) = coolH2Oixd1(:)
      ixd2(:) = coolH2Oixd2(:)
      ixd3(:) = coolH2Oixd3(:)

      !local variables
      v3 = num2col(n(idx_H2O),n(:)) !H2O column density
      cH = n(idx_H) + n(idx_H2)

      v2 = cH
      v1 = inTgas !Tgas

      !logs of variables
      v1 = log10(v1)
      v2 = log10(v2)
      v3 = log10(v3)

      !default value erg/s/cm3
      cooling_H2O = 0d0

      !check limits
      if(v1 >= v1max) v1 = v1max * (1d0 - eps)
      if(v2 >= v2max) v2 = v2max * (1d0 - eps)
      if(v3 >= v3max) v3 = v3max * (1d0 - eps)

      if(v1 < v1min) return
      if(v2 < v2min) return
      if(v3 < v3min) return

      !gets position of variable in the array
      i = (v1 - v1min) * coolH2Odvn1 + 1
      j = (v2 - v2min) * coolH2Odvn2 + 1
      k = (v3 - v3min) * coolH2Odvn3 + 1

      !precompute shared variables
      prev1 = (v1 - x1(i)) * ixd1(i)
      prev2 = (v2 - x2(j)) * ixd2(j)

      !linear interpolation on x1 for x2,x3
      vv1 = prev1 * (coolH2Oy(k, j, i+1) - &
           coolH2Oy(k,j,i)) + coolH2Oy(k,j,i)
      !linear interpolation on x1 for x2+dx2,x3
      vv2 = prev1 * (coolH2Oy(k, j+1, i+1) - &
           coolH2Oy(k, j+1, i)) + coolH2Oy(k, j+1, i)
      !linear interpolation on x2 for x3
      vv12 = prev2 * (vv2 - vv1) + vv1

      !linear interpolation on x1 for x2,x3+dx3
      vv3 = prev1 * (coolH2Oy(k+1, j, i+1) - &
           coolH2Oy(k+1, j, i)) + coolH2Oy(k+1, j, i)
      !linear interpolation on x1 for x2+dx2,x3+dx3
      vv4 = prev1 * (coolH2Oy(k+1, j+1, i+1) - &
           coolH2Oy(k+1,j+1,i)) + coolH2Oy(k+1, j+1, i)
      !linear interpolation on x2 for x3+dx3
      vv34 = prev2 * (vv4 - vv3) + vv3

      !linear interpolation on x3
      xLd = (v3 - x3(k)) * ixd3(k) * (vv34 - &
           vv12) + vv12

      !H2O cooling in erg/s/cm3
      cooling_H2O = 1d1**xLd * cH * n(idx_H2O)

    end function cooling_H2O

    ! ************************
    subroutine init_coolingH2O()
      use krome_commons
      implicit none
      integer::ios,iout(3),i,unit
      real*8::rout(4)

      if (krome_mpi_rank<=1) print *,"load H2O cooling..."
      open(newunit=unit, file="coolH2O.dat", status="old", iostat=ios)
      ! check if file exists
      if(ios.ne.0) then
         print *, "ERROR: problems loading coolH2O.dat!"
         stop
      end if

      do
         read(unit,*,iostat=ios) iout(:),rout(:) !read line
         if(ios<0) exit !eof
         if(ios/=0) cycle !skip blanks
         coolH2Ox1(iout(1)) = rout(1)
         coolH2Ox2(iout(2)) = rout(2)
         coolH2Ox3(iout(3)) = rout(3)
         coolH2Oy(iout(3), iout(2), iout(1)) = rout(4)
      end do

      !store inverse of the differences
      ! to speed up interpolation
      do i=1,coolH2On1-1
         coolH2Oixd1(i) = 1d0 / (coolH2Ox1(i+1) - coolH2Ox1(i))
      end do
      do i=1,coolH2On2-1
         coolH2Oixd2(i) = 1d0 / (coolH2Ox2(i+1) - coolH2Ox2(i))
      end do
      do i=1,coolH2On3-1
         coolH2Oixd3(i) = 1d0 / (coolH2Ox3(i+1) - coolH2Ox3(i))
      end do

      coolH2Ox1min = minval(coolH2Ox1)
      coolH2Ox1max = maxval(coolH2Ox1)
      coolH2Ox2min = minval(coolH2Ox2)
      coolH2Ox2max = maxval(coolH2Ox2)
      coolH2Ox3min = minval(coolH2Ox3)
      coolH2Ox3max = maxval(coolH2Ox3)

      coolH2Odvn1 = (coolH2On1 - 1) / (coolH2Ox1max - coolH2Ox1min)
      coolH2Odvn2 = (coolH2On2 - 1) / (coolH2Ox2max - coolH2Ox2min)
      coolH2Odvn3 = (coolH2On3 - 1) / (coolH2Ox3max - coolH2Ox3min)

    end subroutine init_coolingH2O
#ENDIFKROME

#IFKROME_useCoolingHCN
    !***************************
    ! HCN cooling: courtesy of K.Omukai (Mar2019)
    ! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
    ! see also Omukai+2010 (bit.ly/1HIaGcn)
    ! H and H2 collisions
    function cooling_HCN(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      implicit none
      integer,parameter::imax=coolHCNn1
      integer,parameter::jmax=coolHCNn2
      integer,parameter::kmax=coolHCNn3
      integer::i,j,k
      real*8,parameter::eps=1d-5
      real*8::cooling_HCN,n(:),inTgas
      real*8::v1,v2,v3,prev1,prev2,cH
      real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
      real*8::x1(imax),x2(jmax),x3(kmax)
      real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
      real*8::v1min,v1max,v2min,v2max,v3min,v3max

      !local copy of limits
      v1min = coolHCNx1min
      v1max = coolHCNx1max
      v2min = coolHCNx2min
      v2max = coolHCNx2max
      v3min = coolHCNx3min
      v3max = coolHCNx3max

      !local copy of variables arrays
      x1(:) = coolHCNx1(:)
      x2(:) = coolHCNx2(:)
      x3(:) = coolHCNx3(:)

      ixd1(:) = coolHCNixd1(:)
      ixd2(:) = coolHCNixd2(:)
      ixd3(:) = coolHCNixd3(:)

      !local variables
      v3 = num2col(n(idx_HCN),n(:)) !HCN column density
      cH = n(idx_H) + n(idx_H2)

      v2 = cH
      v1 = inTgas !Tgas

      !logs of variables
      v1 = log10(v1)
      v2 = log10(v2)
      v3 = log10(v3)

      !default value erg/s/cm3
      cooling_HCN = 0d0

      !check limits
      if(v1 >= v1max) v1 = v1max * (1d0 - eps)
      if(v2 >= v2max) v2 = v2max * (1d0 - eps)
      if(v3 >= v3max) v3 = v3max * (1d0 - eps)

      if(v1 < v1min) return
      if(v2 < v2min) return
      if(v3 < v3min) return

      !gets position of variable in the array
      i = (v1 - v1min) * coolHCNdvn1 + 1
      j = (v2 - v2min) * coolHCNdvn2 + 1
      k = (v3 - v3min) * coolHCNdvn3 + 1

      !precompute shared variables
      prev1 = (v1 - x1(i)) * ixd1(i)
      prev2 = (v2 - x2(j)) * ixd2(j)

      !linear interpolation on x1 for x2,x3
      vv1 = prev1 * (coolHCNy(k, j, i+1) - &
           coolHCNy(k,j,i)) + coolHCNy(k,j,i)
      !linear interpolation on x1 for x2+dx2,x3
      vv2 = prev1 * (coolHCNy(k, j+1, i+1) - &
           coolHCNy(k, j+1, i)) + coolHCNy(k, j+1, i)
      !linear interpolation on x2 for x3
      vv12 = prev2 * (vv2 - vv1) + vv1

      !linear interpolation on x1 for x2,x3+dx3
      vv3 = prev1 * (coolHCNy(k+1, j, i+1) - &
           coolHCNy(k+1, j, i)) + coolHCNy(k+1, j, i)
      !linear interpolation on x1 for x2+dx2,x3+dx3
      vv4 = prev1 * (coolHCNy(k+1, j+1, i+1) - &
           coolHCNy(k+1,j+1,i)) + coolHCNy(k+1, j+1, i)
      !linear interpolation on x2 for x3+dx3
      vv34 = prev2 * (vv4 - vv3) + vv3

      !linear interpolation on x3
      xLd = (v3 - x3(k)) * ixd3(k) * (vv34 - &
           vv12) + vv12

      !HCN cooling in erg/s/cm3
      cooling_HCN = 1d1**xLd * cH * n(idx_HCN)

    end function cooling_HCN

    ! ************************
    subroutine init_coolingHCN()
      use krome_commons
      implicit none
      integer::ios,iout(3),i,unit
      real*8::rout(4)

      if (krome_mpi_rank<=1) print *,"load HCN cooling..."
      open(newunit=unit, file="coolHCN.dat", status="old", iostat=ios)
      ! check if file exists
      if(ios.ne.0) then
         print *, "ERROR: problems loading coolHCN.dat!"
         stop
      end if

      do
         read(unit,*,iostat=ios) iout(:),rout(:) !read line
         if(ios<0) exit !eof
         if(ios/=0) cycle !skip blanks
         coolHCNx1(iout(1)) = rout(1)
         coolHCNx2(iout(2)) = rout(2)
         coolHCNx3(iout(3)) = rout(3)
         coolHCNy(iout(3), iout(2), iout(1)) = rout(4)
      end do

      !store inverse of the differences
      ! to speed up interpolation
      do i=1,coolHCNn1-1
         coolHCNixd1(i) = 1d0 / (coolHCNx1(i+1) - coolHCNx1(i))
      end do
      do i=1,coolHCNn2-1
         coolHCNixd2(i) = 1d0 / (coolHCNx2(i+1) - coolHCNx2(i))
      end do
      do i=1,coolHCNn3-1
         coolHCNixd3(i) = 1d0 / (coolHCNx3(i+1) - coolHCNx3(i))
      end do

      coolHCNx1min = minval(coolHCNx1)
      coolHCNx1max = maxval(coolHCNx1)
      coolHCNx2min = minval(coolHCNx2)
      coolHCNx2max = maxval(coolHCNx2)
      coolHCNx3min = minval(coolHCNx3)
      coolHCNx3max = maxval(coolHCNx3)

      coolHCNdvn1 = (coolHCNn1 - 1) / (coolHCNx1max - coolHCNx1min)
      coolHCNdvn2 = (coolHCNn2 - 1) / (coolHCNx2max - coolHCNx2min)
      coolHCNdvn3 = (coolHCNn3 - 1) / (coolHCNx3max - coolHCNx3min)

    end subroutine init_coolingHCN
#ENDIFKROME

#IFKROME_useCoolingZCIENOUV
    !***************************
    !Metal line cooling CIE
    ! method: CLOUDY 10, NOUV
    ! tables kindly provided by Sijing Shen.
    function cooling_Z_CIENOUV(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_fit
      use krome_getphys
      implicit none
      real*8::cooling_Z_CIENOUV,n(:),inTgas
      real*8::cH,Tgas,xLd,logcH

      cooling_Z_CIENOUV = 0d0
      cH = get_Hnuclei(n(:))

      !check if the abundance is close to zero to
      !avoid weird log evaluation
      if(cH.lt.1d-20)return

      Tgas = log10(inTgas)
      logcH = log10(cH)

      xLd = fit_anytab2D(CoolZNOUV_x(:), CoolZNOUV_y(:), CoolZNOUV_z(:,:), &
           CoolZNOUV_xmul, CoolZNOUV_ymul,logcH,Tgas)

      cooling_Z_CIENOUV = 10**xLd * cH * cH * total_Z

    end function cooling_Z_CIENOUV
#ENDIFKROME

#IFKROME_useCoolingZCIE_function
    !***************************
    !Metal line cooling CIE
    ! method: CLOUDY 10, including the
    ! extragalactic (quasars + galaxies) UV flux
    ! by Haardt&Madau 2012.
    !Tables kindly provided by Sijing Shen.
    function cooling_Z_CIE(n,inTgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      implicit none
      integer,parameter::imax=coolZCIEn1
      integer,parameter::jmax=coolZCIEn2
      integer,parameter::kmax=coolZCIEn3
      integer::i,j,k
      real*8::cooling_Z_CIE,n(:),inTgas,Tgas
      real*8::v1,v2,v3,prev1,prev2,cH
      real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
      real*8::x1(imax),x2(jmax),x3(kmax)
      real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
      real*8::v1min,v1max,v2min,v2max,v3min,v3max
      real*8,parameter::eps=1d-5

      Tgas = inTgas
      cooling_Z_CIE = 0d0

      !local copy of limits
      v1min = coolZCIEx1min
      v1max = coolZCIEx1max
      v2min = coolZCIEx2min
      v2max = coolZCIEx2max
      v3min = coolZCIEx3min
      v3max = coolZCIEx3max

      !local copy of variables arrays
      x1(:) = coolZCIEx1(:)
      x2(:) = coolZCIEx2(:)
      x3(:) = coolZCIEx3(:)

      ixd1(:) = coolZCIEixd1(:)
      ixd2(:) = coolZCIEixd2(:)
      ixd3(:) = coolZCIEixd3(:)

      !local variables
      cH = get_Hnuclei(n(:))

      !check if the abundance is close to zero to
      !avoid weird log evaluation
      if(cH.lt.1d-20)return

      v1 = Tgas           !Tgas
      v2 = cH             !total H number density
      v3 = phys_zredshift !redshift is linear

      !logs of variables
      v1 = log10(v1)
      v2 = log10(v2)

      !check limits
      if(v1>=v1max) v1 = v1max*(1d0-eps)
      if(v2>=v2max) v2 = v2max*(1d0-eps)
      if(v3>=v3max) v3 = v3max*(1d0-eps)

      if(v1<v1min) return
      if(v2<v2min) return
      if(v3<v3min) return

      !gets position of variable in the array
      i = (v1-v1min)*coolZCIEdvn1+1
      j = (v2-v2min)*coolZCIEdvn2+1
      k = (v3-v3min)*coolZCIEdvn3+1

      !precompute shared variables
      prev1 = (v1-x1(i))*ixd1(i)
      prev2 = (v2-x2(j))*ixd2(j)

      !linear interpolation on x1 for x2,x3
      vv1 = prev1 * (coolZCIEy(k,j,i+1) - &
           coolZCIEy(k,j,i)) + coolZCIEy(k,j,i)
      !linear interpolation on x1 for x2+dx2,x3
      vv2 = prev1 * (coolZCIEy(k,j+1,i+1) - &
           coolZCIEy(k,j+1,i)) + coolZCIEy(k,j+1,i)
      !linear interpolation on x2 for x3
      vv12 = prev2 * (vv2 - vv1) + vv1

      !linear interpolation on x1 for x2,x3+dx3
      vv3 = prev1 * (coolZCIEy(k+1,j,i+1) - &
           coolZCIEy(k+1,j,i)) + coolZCIEy(k+1,j,i)
      !linear interpolation on x1 for x2+dx2,x3+dx3
      vv4 = prev1 * (coolZCIEy(k+1,j+1,i+1) - &
           coolZCIEy(k+1,j+1,i)) + coolZCIEy(k+1,j+1,i)
      !linear interpolation on x2 for x3+dx3
      vv34 = prev2 * (vv4 - vv3) + vv3

      !linear interpolation on x3
      xLd = (v3-x3(k))*ixd3(k)*(vv34 - &
           vv12) + vv12

      !Z cooling in erg/s/cm3
      cooling_Z_CIE = 1d1**xLd * cH * cH * total_Z

    end function cooling_Z_CIE

    !************************
    subroutine init_coolingZCIE()
      use krome_commons
      implicit none
      integer::ios,iout(3),i
      real*8::rout(5)

      if(krome_mpi_rank<=1) print *,"load Z_CIE2012 cooling..."
      open(33,file="coolZ_CIE2012.dat",status="old",iostat=ios)
      !check if file exists
      if(ios.ne.0) then
         print *,"ERROR: problems loading coolZ_CIE2012.dat!"
         stop
      end if

      do
         read(33,*,iostat=ios) iout(:),rout(:) !read line
         if(ios<0) exit !eof
         if(ios/=0) cycle !skip blanks
         coolZCIEx1(iout(1)) = rout(1)
         coolZCIEx2(iout(2)) = rout(2)
         coolZCIEx3(iout(3)) = rout(3)
         coolZCIEy(iout(3),iout(2),iout(1)) = rout(4)
         heatZCIEy(iout(3),iout(2),iout(1)) = rout(5)
      end do

      !store inverse of the differences
      ! to speed up interpolation
      do i=1,coolZCIEn1-1
         coolZCIEixd1(i) = 1d0/(coolZCIEx1(i+1)-coolZCIEx1(i))
      end do
      do i=1,coolZCIEn2-1
         coolZCIEixd2(i) = 1d0/(coolZCIEx2(i+1)-coolZCIEx2(i))
      end do
      do i=1,coolZCIEn3-1
         coolZCIEixd3(i) = 1d0/(coolZCIEx3(i+1)-coolZCIEx3(i))
      end do

      coolZCIEx1min = minval(coolZCIEx1)
      coolZCIEx1max = maxval(coolZCIEx1)
      coolZCIEx2min = minval(coolZCIEx2)
      coolZCIEx2max = maxval(coolZCIEx2)
      coolZCIEx3min = minval(coolZCIEx3)
      coolZCIEx3max = maxval(coolZCIEx3)

      coolZCIEdvn1 = (coolZCIEn1-1)/(coolZCIEx1max-coolZCIEx1min)
      coolZCIEdvn2 = (coolZCIEn2-1)/(coolZCIEx2max-coolZCIEx2min)
      coolZCIEdvn3 = (coolZCIEn3-1)/(coolZCIEx3max-coolZCIEx3min)

    end subroutine init_coolingZCIE
#ENDIFKROME

#IFKROME_useCoolingZExtended
    !************************
    function cooling_ZExtended(n,Tgas)
      use krome_commons
      implicit none
      real*8::cooling_ZExtended
      real*8::n(:),Tgas
      real*8::f1,f2,smooth

      !this parameter controls the smoothness of the
      ! merge between the two cooling functions
      smooth = 1.d-3

      !smoothing functions | f1+f2=1
      f1 = (tanh(smooth*(Tgas-1d4))+1.d0)*0.5d0
      f2 = (tanh(smooth*(-Tgas+1d4))+1.d0)*0.5d0

      cooling_ZExtended = f1*cooling_Z_CIE(n(:),Tgas) &
           + f2*(cooling_Z(n(:),Tgas) #KROME_floorZ_EXTendED)

    end function cooling_ZExtended
#ENDIFKROME


    !*****************************
    function cooling_custom(n,Tgas)
      use krome_commons
      use krome_subs
      use krome_constants
      implicit none
      real*8::n(:),Tgas,cooling_custom
#KROME_custom_cooling_var_define

      cooling_custom = 0d0
#KROME_custom_cooling_var
#KROME_custom_cooling_expr

    end function cooling_custom

    !**********************************
    function kpla(n,Tgas)
      !Planck opacity mean fit (Lenzuni+1996)
      !only denisity dependent (note that the
      ! fit provided by Lenzuni is wrong)
      ! valid for T<3e3 K
      !use krome_subs
      use krome_commons
      use krome_getphys
      implicit none
      real*8::kpla,rhogas,Tgas,n(:),y
      real*8::a0,a1,m(nspec)

      m(:) = get_mass()
      rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

      kpla = 0.d0
      !opacity is zero under 1e-12 g/cm3
      if(rhogas<1d-12) return

      !fit coefficients
      a0 = 1.000042d0
      a1 = 2.14989d0

      !log density cannot exceed 0.5 g/cm3
      y = log10(min(rhogas,0.5d0))

      kpla = 1d1**(a0*y + a1) !fit density only

    end function kpla

    !*****************************
    function coolingChem(n,Tgas)
      implicit none
      real*8::coolingChem,n(:),Tgas

      !note that this function is a dummy.
      ! For chemical cooling you should see
      ! heatingChem function in krome_heating.f90

      coolingChem = 0.d0

    end function coolingChem

#IFKROME_useCoolingContinuum
    !**********************************
    function cooling_Continuum(n,Tgas)
      !cooling from continuum for a thin gas (no opacity)
      !see Omukai+2000 for details
      use krome_commons
      use krome_constants
      use krome_subs
      use krome_getphys
      implicit none
      real*8::n(:),Tgas,cooling_Continuum,kgas,rhogas
      real*8::lj,tau,beta,m(nspec)

      m(:) = get_mass()
      rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3
      kgas = kpla(n(:),Tgas) !planck opacity cm2/g (Omukai+2000)
      lj = get_jeans_length(n(:), Tgas) !cm
      tau = lj * kgas * rhogas + 1d-40 !opacity
      beta = min(1.d0,tau**(-2)) !beta escape (always <1.)
      cooling_Continuum = 4.d0 * stefboltz_erg * Tgas**4 &
           * kgas * rhogas * beta !erg/s/cm3

    end function cooling_Continuum
#ENDIFKROME


#IFKROME_useCoolingCIE
    !*******************************
    function cooling_CIE(n, inTgas)
      !CIE cooling: fit from Ripamponti&Abel2004 (RA04) data
      ! The fit is valid from 100K-1e6K.
      ! Original data are from 400K to 7000K.
      ! We extrapolated data under 400K and fitted from 100K to 10**2.95 K.
      ! Data from 10**2.95 K to 1e5K are fitted analogously.
      ! Above 1e5 we employ a cubic extrapolation.
      use krome_commons
      use krome_constants
      real*8::cooling_CIE,n(:),Tgas,inTgas
      real*8::x,x2,x3,x4,x5
      real*8::a0,a1,a2,a3,a4,a5
      real*8::b0,b1,b2,b3,b4,b5
      real*8::cool,tauCIE,logcool

      cooling_CIE = 0d0
      !set cooling to zero if n_H2 is smaller than 1e-12 1/cm3
      !to avoid division by zero in opacity term due to tauCIE=0
      if(n(idx_H2)<1d-12) return

      Tgas = inTgas
      !temperature limit
      if(Tgas<phys_Tcmb) return

      !prepares variables
      x = log10(Tgas)
      x2 = x*x
      x3 = x2*x
      x4 = x3*x
      x5 = x4*x

      cool = 0.d0
      !outside boundaries below cooling is zero
      logcool = -1d99

      !evaluates fitting functions
      if(x>2.d0 .and. x<2.95d0) then
         a0 = -30.3314216559651d0
         a1 = 19.0004016698518d0
         a2 = -17.1507937874082d0
         a3 = 9.49499574218739d0
         a4 = -2.54768404538229d0
         a5 = 0.265382965410969d0
         logcool = a0 + a1*x + a2*x2 + a3*x3 +a4*x4 +a5*x5
      elseif(x.GE.2.95d0 .and. x<5.d0) then
         b0 = -180.992524120965d0
         b1 = 168.471004362887d0
         b2 = -67.499549702687d0
         b3 = 13.5075841245848d0
         b4 = -1.31983368963974d0
         b5 = 0.0500087685129987d0
         logcool = b0 + b1*x + b2*x2 + b3*x3 +b4*x4 +b5*x5
      elseif(x.GE.5.d0) then
         logcool = 3.d0 * x - 21.2968837223113 !cubic extrapolation
      end if

      !opacity according to RA04
      tauCIE = (n(idx_H2) * 1.4285714d-16)**2.8 !note: 1/7d15 = 1.4285714d-16
      cool = p_mass * 1d1**logcool !erg*cm3/s

      cooling_CIE = cool * min(1.d0, (1.d0-exp(-tauCIE))/tauCIE) &
           * n(idx_H2) * sum(n(1:nmols)) !erg/cm3/s

    end function cooling_CIE
#ENDIFKROME


#IFKROME_useCoolingExpansion
    !*******************************
    function cooling_expansion(n, Tgas)
      !R'/R expansion cooling erg/cm3/s from Galli&Palla 1998
      use krome_user_commons
      use krome_commons
      use krome_constants
      real*8::cooling_expansion,n(:),Tgas
      real*8::ntot

      ntot=sum(n(1:nmols))
      !note that user_* must be provided by the user
      ! in the reaction file using @common: and initialized
      ! using the interface subroutine
      cooling_expansion = 3.d0*ntot*boltzmann_erg*Tgas*hubble0 &
           * (1.d0 + phys_zredshift) &
           * sqrt(omega0 * phys_zredshift + 1.d0) !erg/s/cm3

    end function cooling_expansion
#ENDIFKROME

#IFKROME_useCoolingCompton
    !*******************************
    function cooling_compton(n, Tgas)
      !compton cooling erg/cm3/s from Cen1992
      use krome_user_commons
      use krome_commons
      real*8::cooling_compton,n(:),Tgas

      !note that redhsift is a common variable and
      ! should be provided by the user, otherwise the default is zero
      cooling_compton = 5.65d-36 * (1.d0 + phys_zredshift)**4 &
           * (Tgas - 2.73d0 * (1.d0 + phys_zredshift)) * n(idx_e) !erg/s/cm3

    end function cooling_compton
#ENDIFKROME

#IFKROME_useCoolingDustTabs
    !****************************
    function cooling_dust(n,Tgas)
      use krome_commons
      use krome_subs
      use krome_fit
      use krome_getphys
      implicit none
      real*8::n(:),Tgas,ntot,cooling_dust,coolFit
      real*8::logn,logt

      ntot = sum(n(1:nmols))
      Tgas = n(idx_Tgas)

      logn = log10(ntot)
      logt = log10(Tgas)
      !cooling fit from tables
#IFKROME_dust_table_2D
      coolFit = fit_anytab2D_linlog(dust_tab_ngas(:), dust_tab_Tgas(:), &
           dust_tab_cool(:,:), dust_mult_ngas, dust_mult_Tgas, &
           logn, logt)
#ENDIFKROME_dust_table_2D

#IFKROME_dust_table_3D
      coolFit = fit_anytab3D_linlinlog(dust_tab_ngas(:), dust_tab_Tgas(:), &
           dust_tab_AvVariable(:), &
           dust_tab_cool(:,:,:), dust_mult_ngas, dust_mult_Tgas, &
           dust_mult_AvVariable, &
           logn, logt, dust_table_AvVariable_log)
#ENDIFKROME_dust_table_3D

      cooling_dust = get_mu(n) * coolFit * ntot * ntot

    end function cooling_dust
#ENDIFKROME

#IFKROME_useCoolingDustNoTdust
    !*******************************
    function cooling_dust(n,Tgas)
      !cooling from dust in erg/cm3/s
      use krome_constants
      use krome_commons
      use krome_dust
      use krome_getphys
      implicit none
      real*8::cooling_dust,n(:),Tgas
      real*8::pre,ntot,vgas,fact
      integer::i
      fact = 0.5d0
      Tgas = n(idx_Tgas)
      vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas
      ntot = sum(n(1:nmols))
      pre = 2d0*fact*vgas*boltzmann_erg*ntot

      cooling_dust = 0d0
      do i=1,ndust
         cooling_dust = cooling_dust &
              + pre * xdust(i) * krome_dust_asize2(i) &
              * (Tgas-krome_dust_T(i))
      end do

    end function cooling_dust
#ENDIFKROME


#IFKROME_useCoolingDust
    !*******************************
    function cooling_dust(n,Tgas)
      !cooling from dust in erg/cm3/s
      use krome_constants
      use krome_commons
      use krome_subs
      use krome_dust
      use krome_getphys
      implicit none
      real*8::cooling_dust,n(:),Tgas
#IFKROME_usedTdust
      real*8::rhogas,ljeans,be,ntot,vgas
      real*8::m(nspec),intCMB,fact,intJflux
      integer::i

      fact = 0.5d0
      cooling_dust = 0d0
      m(:) = get_mass()
      Tgas = n(idx_Tgas)
      vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas
      ntot = sum(n(1:nmols))
      rhogas = sum(n(1:nmols)*m(1:nmols))
      ljeans = get_jeans_length_rho(n(:),Tgas,rhogas)
      be = besc(n(:),Tgas,ljeans,rhogas)
      intJflux = 0d0

      do i=1,ndust
#IFKROME_usePhotoDust
         !compute external radiation term
         intJflux = get_int_JQabs(i)
#ENDIFKROME_usePhotoDust
         intCMB = get_dust_intBB(i,phys_Tcmb)
         cooling_dust = cooling_dust + (get_dust_intBB(i,n(nmols+ndust+i)) &
              - intCMB - intJflux) * be * xdust(i) * krome_dust_asize2(i)
      end do
      cooling_dust = 4d0*pi*cooling_dust !erg/s/cm3
      return
#ENDIFKROME_usedTdust

      cooling_dust = dust_cooling !erg/s/cm3

    end function cooling_dust
#ENDIFKROME


#IFKROME_useCoolingdH
    !*******************************
    function cooling_dH(n,Tgas)
      !cooling from reaction enthalpy erg/s/cm3
      use krome_commons
      implicit none
      real*8::cooling_dH,cool,n(:),Tgas,small,nmax
      real*8::logT,lnT,Te,lnTe,T32,t3,invT,invTe,sqrTgas,invsqrT32,sqrT32
#KROME_vars

      !replace small according to the desired enviroment
      ! and remove nmax if needed
      nmax = maxval(n(1:nmols))
      small = #KROME_small

      logT = log10(Tgas) !log10 of Tgas (#)
      lnT = log(Tgas) !ln of Tgas (#)
      Te = Tgas*8.617343d-5 !Tgas in eV (eV)
      lnTe = log(Te) !ln of Te (#)
      T32 = Tgas/3.d2 !Tgas/(300 K) (#)
      t3 = T32 !alias for T32 (#)
      invT = 1.d0/Tgas !inverse of T (1/K)
      invTe = 1.d0/Te !inverse of T (1/eV)
      sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)
      invsqrT32 = 1.d0/sqrt(T32)
      sqrT32 = sqrt(T32)

      cool = 0.d0

#KROME_rates
#KROME_dH_cooling

      cooling_dH = cool

    end function cooling_dH
#ENDIFKROME

#IFKROME_useH2esc_omukai
    !*****************************
    !escape opacity for H2 cooling.
    ! courtesy of Kazu Omukai (2014)
    ! Einstein's A coefficients for spontaneous emission
    ! calculated by Turner, Kirby-Docken, & Dalgarno 1977, ApJS, 35, 281
    ! and the excitation energies for the levels of Borysow,
    ! Frommhold & Moraldi (1989), ApJ, 336, 495.
    function H2opacity_omukai(Tgas, n)
      use krome_commons
      use krome_fit
      use krome_subs
      use krome_getphys
      implicit none
      real*8::H2opacity_omukai,Tgas,ntot,lTgas,lntot,n(:)

      ntot = sum(n(1:nmols))
      lTgas = log10(Tgas)
      lntot = log10(num2col(ntot,n(:)))

      H2opacity_omukai = 1d1**(fit_anytab2D(arrH2esc_ntot(:), &
           arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
           ymulH2esc,lntot,lTgas))

    end function H2opacity_omukai
#ENDIFKROME

#IFKROME_useCoolingH2GP
    !*******************************
    function cooling_H2GP(n, Tgas)
      !cooling from Galli&Palla98
      use krome_commons
      use krome_subs
      real*8::n(:),Tgas, tm, logT
      real*8::cooling_H2GP,T3
      real*8::LDL,HDLR,HDLV,HDL

      tm = max(Tgas, 13.0d0)    ! no cooling below 13 Kelvin
      tm = min(Tgas, 1.d5)      ! fixes numerics
      logT = log10(tm)
      T3 = tm * 1.d-3

      !low density limit in erg/s
      LDL = 1.d1**(-103.d0+97.59d0*logT-48.05d0*logT**2&
           +10.8d0*logT**3-0.9032d0*logT**4)*n(idx_H)

      !this will avoid a division by zero and useless calculations
      if(LDL==0d0) then
         cooling_H2GP = 0d0
         return
      end if

      !high density limit
      HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
           3.e-24*exp(-0.51/t3)) !erg/s
      HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
      HDL  = HDLR + HDLV !erg/s

      !to avoid division by zero
      if (HDL==0d0) then
         cooling_H2GP = 0d0
      else
         cooling_H2GP = n(idx_H2)/(1d0/HDL+1d0/LDL) #KROME_H2opacity !erg/cm3/s
      endif

    end function cooling_H2GP
#ENDIFKROME

#IFKROME_useCoolingH2

    !*****************
    !sigmoid function with x0 shift and s steepness
    function sigmoid(x,x0,s)
      implicit none
      real*8::sigmoid,x,x0,s

      sigmoid = 1d1/(1d1+exp(-s*(x-x0)))

    end function sigmoid

    !*******************
    !window function for H2 cooling to smooth limits
    function wCool(logTgas,logTmin,logTmax)
      implicit none
      real*8::wCool,logTgas,logTmin,logTmax,x

      x = (logTgas-logTmin)/(logTmax-logTmin)
      wCool = 1d1**(2d2*(sigmoid(x,-2d-1,5d1)*sigmoid(-x,-1.2d0,5d1)-1d0))
      if(wCool<1d-199) wCool = 0d0
      if(wCool>1d0) then
         print *,"ERROR: wCool>1"
         stop
      end if

    end function wCool

    !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
    !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
    !UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
    !EACH SINGLE FUNCTION IS IN erg/s
    !FINAL UNITS = erg/cm3/s
    !*******************************
    function cooling_H2(n, Tgas)
      use krome_commons
      use krome_subs
      use krome_getphys
      real*8::n(:),Tgas
      real*8::temp,logt3,logt,cool,cooling_H2,T3
      real*8::LDL,HDLR,HDLV,HDL
      real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
      real*8::dump14,fH2H,fH2e,fH2H2,fH2Hp,fH2He,w14,w24
      integer::i
      character*16::names(nspec)

      temp = Tgas
      cooling_H2 = 0d0
      !if(temp<2d0) return

      T3 = temp * 1.d-3
      logt3 = log10(T3)
      logt = log10(temp)
      cool = 0d0

      logt32 = logt3 * logt3
      logt33 = logt32 * logt3
      logt34 = logt33 * logt3
      logt35 = logt34 * logt3
      logt36 = logt35 * logt3
      logt37 = logt36 * logt3
      logt38 = logt37 * logt3

      w14 = wCool(logt, 1d0, 4d0)
      w24 = wCool(logt, 2d0, 4d0)

#IFKROME_hasH
      !//H2-H
      if(temp<=1d2) then
         fH2H = 1.d1**(-16.818342D0 +3.7383713D1*logt3 &
              +5.8145166D1*logt32 +4.8656103D1*logt33 &
              +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
      elseif(temp>1d2 .and. temp<=1d3) then
         fH2H = 1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
              -1.1332860D1*logt32 -2.7850082D1*logt33 &
              -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
      elseif(temp>1.d3.and.temp<=6d3) then
         fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
              -3.7209846D0*logt32 +5.9369081D0*logt33 &
              -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
      else
         fH2H = 1.862314467912518E-022*wCool(logt,1d0,log10(6d3))*n(idx_H)
      end if
      cool = cool + fH2H
#ENDIFKROME_hasH

#IFKROME_hasHp
      !//H2-Hp
      if(temp>1d1.and.temp<=1d4) then
         fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
              +0.015391166d0*logt32 -0.23619985d0*logt33 &
              -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj)
      else
         fH2Hp = 1.182509139382060E-021*n(idx_Hj)*w14
      endif
      cool = cool + fH2Hp
#ENDIFKROME_hasHp

      !//H2-H2
      fH2H2 = w24*1d1**(-2.3962112D1 +2.09433740D0*logt3 &
           -.77151436D0*logt32 +.43693353D0*logt33 &
           -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) !&
      cool = cool + fH2H2

#IFKROME_hasElectrons
      !//H2-e
      fH2e = 0d0
      if(temp<=5d2) then
         fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
              +9.6743155d1*logt32 +3.4319180d2*logt33 &
              +7.3471651d2*logt34 +9.8367576d2*logt35 &
              +8.0181247d2*logt36 +3.6414446d2*logt37 &
              +7.0609154d1*logt38,3d1))*n(idx_e)
      elseif(temp>5d2)  then
         fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
              +.93310622D0*logt32 +4.0406627d0*logt33 &
              -4.7274036d0*logt34 -8.8077017d0*logt35 &
              +8.9167183*logt36 + 6.4380698*logt37 &
              -6.3701156*logt38)*n(idx_e)
      end if
      cool = cool + fH2e*w24
#ENDIFKROME_hasElectrons

#IFKROME_hasHe
      !//H2-He
      if(temp>1d1.and.temp<=1d4)then
         fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3&
              -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
              +.19191375d0*logt35)*n(idx_He)
      else
         fH2He = 1.002560385050777E-022*n(idx_He)*w14
      endif
      cool = cool + fH2He
#ENDIFKROME_hasHe

      !check error
      if(cool>1.d30) then
         print *,"ERROR: cooling >1.d30 erg/s/cm3"
         print *,"cool (erg/s/cm3): ",cool
         names(:) = get_names()
         do i=1,size(n)
            print '(I3,a18,E11.3)',i,names(i),n(i)
         end do
         stop
      end if

      !this to avoid negative, overflow and useless calculations below
      if(cool<=0d0) then
         cooling_H2 = 0d0
         return
      end if

      !high density limit from HM79, GP98 below Tgas = 2d3
      !UPDATED USING GLOVER 2015 for high temperature corrections, MNRAS
      !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
      !the following mix of functions ensures the right behaviour
      ! at low (T<10 K) and high temperatures (T>2000 K) by
      ! using both the original Hollenbach and the new Glover data
      ! merged in a smooth way.
      if(temp.lt.2d3)then
         HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
              3.e-24*exp(-0.51/t3)) !erg/s
         HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
         HDL  = HDLR + HDLV !erg/s
      elseif(temp>=2d3 .and. temp<=1d4)then
         HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
              -1.5738805*logt32 -4.7155769*logt33 &
              +2.4714161*logt34 +5.4710750*logt35 &
              -3.9467356*logt36 -2.2148338*logt37 &
              +1.8161874*logt38)
      else
         dump14 = 1d0 / (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
         HDL = 5.531333679406485E-019*dump14
      endif

      LDL = cool !erg/s
      if (HDL==0.) then
         cooling_H2 = 0.d0
      else
         cooling_H2 = n(idx_H2)/(1.d0/HDL+1.d0/LDL) #KROME_H2opacity !erg/cm3/s
      endif

    end function cooling_H2
#ENDIFKROME


#IFKROME_useCoolingAtomic
    !Atomic COOLING  Cen ApJS, 78, 341, 1992
    !UNITS = erg/s/cm3
    !*******************************
    function cooling_Atomic(n, Tgas)
      use krome_commons
      use krome_subs
      real*8::Tgas,cooling_atomic,n(:)
      real*8::temp,T5,cool


      temp = max(Tgas,10d0) !K
      T5 = temp/1d5 !K
      cool = 0d0 !erg/cm3/s

      !COLLISIONAL IONIZATION: H, He, He+, He(2S)
#IFKROME_hasH
      cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
           *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)
#ENDIFKROME_hasH

#IFKROME_hasHe
      cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
           *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)
#ENDIFKROME_hasHe


#IFKROME_hasHep
      cool = cool+ 4.95d-22*sqrt(temp)/(1.d0+sqrt(T5))&
           *exp(-6.31515d5/temp)*n(idx_e)*n(idx_Hej)
      cool = cool+ 5.01d-27*temp**(-0.1687)/(1.d0+sqrt(T5))&
           *exp(-5.5338d4/temp)*n(idx_e)**2*n(idx_Hej)
#ENDIFKROME_hasHep

      !RECOMBINATION: H+, He+,He2+
#IFKROME_hasHp
      cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
           /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)
#ENDIFKROME_hasHp

#IFKROME_hasHep
      cool = cool+ 1.55d-26*temp**(0.3647)*n(idx_e)*n(idx_Hej)
#ENDIFKROME_hasHp

#IFKROME_hasHepp
      cool = cool+ 3.48d-26*sqrt(temp)*(temp/1.d3)**(-0.2)&
           /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hejj)
#ENDIFKROME_hasHepp

#IFKROME_hasHe
      !DIELECTRONIC RECOMBINATION: He
      cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
           *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)
#ENDIFKROME_hasHe

      !COLLISIONAL EXCITATION:
      !H(all n), He(n=2,3,4 triplets), He+(n=2)
#IFKROME_hasH
      cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)
#ENDIFKROME_hasH

#IFKROME_hasHep
      cool = cool+ 9.1d-27*temp**(-.1687)/(1.d0+sqrt(T5))&
           *exp(-1.3179d4/temp)*n(idx_e)**2*n(idx_Hej)
      cool = cool+ 5.54d-17*temp**(-.397)/(1.d0+sqrt(T5))&
           *exp(-4.73638d5/temp)*n(idx_e)*n(idx_Hej)
#ENDIFKROME_hasHep

      cooling_atomic = max(cool, 0d0)  !erg/cm3/s

    end function cooling_Atomic
#ENDIFKROME

#IFKROME_useCoolingFF

    !**************************
    !free-free cooling (bremsstrahlung for all ions)
    ! using mean Gaunt factor value (Cen+1992)
    function cooling_ff(n,Tgas)
      use krome_commons
      implicit none
      real*8::n(:),Tgas,cool,cooling_ff,gaunt_factor,bms_ions

      gaunt_factor = 1.5d0 !mean value

      !BREMSSTRAHLUNG: all ions
#KROME_brem_ions
      cool = 1.42d-27*gaunt_factor*sqrt(Tgas)&
           *bms_ions*n(idx_e)

      cooling_ff = max(cool, 0.d0)  !erg/cm3/s

    end function cooling_ff
#ENDIFKROME

#IFKROME_useCoolingHD
    !HD COOLING LIPOVKA ET AL. MNRAS, 361, 850, (2005)
    !UNITS=erg/cm3/s
    !*******************************
    function cooling_HD(n, inTgas)
      use krome_commons
      use krome_subs
      implicit none
      integer::i,j
      integer, parameter::ns=4
      real*8::cooling_HD
      real*8::n(:),Tgas,logTgas,lognH,inTgas
      real*8::c(0:ns,0:ns),logW,W,dd,lhj

      !default HD cooling value
      cooling_HD = 0.0d0 !erg/cm3/s

      !this function does not have limits on density
      ! and temperature, even if the original paper do.
      ! However, we extrapolate the limits.

      !exit on low temperature
      if(inTgas<phys_Tcmb) return
      !extrapolate higher temperature limit
      Tgas = min(inTgas,1d4)

      !calculate density
      dd = n(idx_H) !sum(n(1:nmols))
      !exit if density is out of Lipovka bounds (uncomment if needed)
      !if(dd<1d0 .or. dd>1d8) return

      !extrapolate density limits
      dd = min(max(dd,1d-2),1d10)

      !POLYNOMIAL COEFFICIENT: TABLE 1 LIPOVKA
      c(0,:) = (/-42.56788d0, 0.92433d0, 0.54962d0, -0.07676d0, 0.00275d0/)
      c(1,:) = (/21.93385d0, 0.77952d0, -1.06447d0, 0.11864d0, -0.00366d0/)
      c(2,:) = (/-10.19097d0, -0.54263d0, 0.62343d0, -0.07366d0, 0.002514d0/)
      c(3,:) = (/2.19906d0, 0.11711d0, -0.13768d0, 0.01759d0, -0.00066631d0/)
      c(4,:) = (/-0.17334d0, -0.00835d0, 0.0106d0, -0.001482d0, 0.00006192d0/)

      logTgas = log10(Tgas)
      lognH   = log10(dd)

      !loop to compute coefficients
      logW = 0.d0
      do j = 0, ns
         lHj = lognH**j
         do i = 0, ns
            logW = logW + c(i,j)*logTgas**i*lHj !erg/s
         enddo
      enddo

      W = 10.d0**(logW)
      cooling_HD = W * n(idx_HD) !erg/cm3/s


    end function cooling_HD
#ENDIFKROME

#IFKROME_useCoolingGH
    !******************************
    function cooling_GH(n, Tgas)
      use krome_commons
      use krome_constants, only : ip_mass
      use krome_coolingGH
      use krome_getphys, only : get_mass, get_metallicityC, get_metallicityO
      implicit none
      real*8::n(:),Tgas
      real*8:: cooling_GH
      real*8:: PLW_last=0d0,PHI_last=0d0,PHeI_last=0d0, &
               PCVI_last=0d0,nb_last=0d0
      real*8,  save :: rch(NRCH)
      integer, save :: ich(NICH)
      !$omp threadprivate(PLW_last,PHI_last,PHeI_last,PCVI_last,nb_last,rch,ich)
      real*8::log10Tgas,nb,Zmetal,ZAsplund,cfun,hfun
      logical, save :: first_call=.true.
      integer       :: ierr
      ! Parameters related to dust destruction at high temperatures
      real*8, parameter :: Tdestruct=1e6, smooth=1./(0.1*Tdestruct)
      ! Check if we have to load table -- only done by one thread
      ! Do the double if-block so that we only enter the critical region
      ! if first_call has not happened or is in progress
      !----------------------------------
      if (first_call) then
         !$omp critical
         if (first_call) then
            ! NOTICE ierr is double duty.
            ! On input: sets mode of tables. =0: full heating nd cooling function. =1: only metals
            ! On output: error code
            ierr = 1
            call frtInitCF(ierr,'cf_table.I2.dat')
            if(ierr .ne. 0) then
               print *,'Error in reading Gnedin and Hollon'
               print *,' cooling table data file: ', ierr
               print *,'Maybe data file cf_table.I2.dat'
               print *,' does not exist in current directory ?'
               stop
            endif
         end if
         first_call = .false.
         !$omp end critical
      end if

      nb = sum(n(1:nspec)*get_mass()) * ip_mass

      if (Tgas > 1d3) then
        !check input values to see if we have to regenerate the cache
        if (abs(nb-nb_last) > 1d-6*nb .or. &
             PLW  .ne. PLW_last  .or. &
             PHI  .ne. PHI_last  .or. &
             PHeI .ne. PHeI_last .or. &
             PCVI .ne. PCVI_last) then
           ! Gas-phase metallicity based on an average of C and O in case they differ wrt dust depletion
           Zasplund = 10.**(0.5*(get_metallicityC(n)+get_metallicityO(n))) ! Normalised to Asplund 2009 value of Z=0.0134
           Zmetal = Zasplund * 0.0134 / 0.02 ! Metallicity of gas in units of 0.02
           ! Above a certain temperature, Tdestruct, correct depletion by
           ! assuming dust-to-gas = 0.01 and all dust sublimated and therefore
           ! Zmetal = (Zmetal * 0.02 + 0.01) / 0.02
           ! At Tdestruct = 10^6 K 1% of protons have v > 150 km/s and can directly destruct dust grains
           Zmetal = (Zmetal * 0.02 + &
                     0.01 * ((tanh(smooth*(Tgas-Tdestruct))+1.d0)*0.5d0) &
                    ) / 0.02
           !
           call frtCFCache(nb, Zmetal, PLW, PHI, PHeI, PCVI, ich, rch, ierr)
           nb_last=nb
           PLW_last=PLW
           PHI_last=PHI
           PHeI_last=PHeI
           PCVI_last=PCVI
           if (ierr > 0) then
              print *,'Problems with caching Gnedin and Hollon cooling table.'
              print *,' Stopping. Error code :', ierr
              print *,'Input variables nb, Tgas, PLW, PHI, PHeI, PCVI: ', nb, Tgas, PLW, PHI, PHeI, PCVI
              stop
           end if
        end if

        call frtCFGetLn(log(max(1d0,Tgas)),ich,rch,cfun,hfun)
        !call frtGetCF(Tgas, nb, 1d0, Plw, Phi, Phei, Pcvi, cfun, hfun, ierr)

      else
        cfun = 0d0
        hfun = 0d0
      endif

      cooling_GH = nb**2 * (cfun-hfun)

    end function cooling_GH
#ENDIFKROME

#IFKROME_useCoolingZ_function
    !*********************************************

    !function for linear interpolation of f(x), using xval(:)
    ! and the corresponding yval(:) as reference values
    ! note: slow function, use only for initializations
    function flin(xval,yval,x)
      implicit none
      real*8::xval(:),yval(:),x,flin
      integer::i,n
      logical::found
      found = .false.
      n = size(xval)
      x = max(x,xval(1)) !set lower bound
      x = min(x,xval(n)) !set upper bound
      !loop to find interval (slow)
      do i=2,n
         if(x.le.xval(i)) then
            !linear fit
            flin = (yval(i) - yval(i-1)) / (xval(i) - xval(i-1)) * &
                 (x - xval(i-1)) + yval(i-1)
            found = .true. !found flag
            exit
         end if
      end do
      if(.not.found) flin = yval(n)

    end function flin

    !************************
    !dump the level populations in a file
    subroutine dump_cooling_pop(Tgas,nfile)
      implicit none
      integer::nfile,i
      real*8::Tgas

#KROME_popvar_dump
      write(nfile,*)

    end subroutine dump_cooling_pop

    !***********************
    !metal cooling as in Maio et al. 2007
    ! loaded from data file
    function cooling_Z(n,inTgas)
      use krome_commons
      use krome_constants
      implicit none
      real*8::n(:), inTgas, cool, cooling_Z, k(nZrate), Tgas

      Tgas = inTgas
      k(:) = coolingZ_rate_tabs(Tgas)

      cool = 0d0
#KROME_coolingZ_call_functions

      cooling_Z = cool * boltzmann_erg

    end function cooling_Z

    !********************************
    function coolingZ_rates(inTgas)
      use krome_commons
      use krome_subs
      use krome_fit
      implicit none
      real*8::inTgas,coolingZ_rates(nZrate),k(nZrate)
      real*8::Tgas,invT,logTgas
      integer::i
#KROME_coolingZ_declare_custom_vars

      Tgas = inTgas
      invT = 1d0/Tgas
      logTgas = log10(Tgas)

#KROME_coolingZ_custom_vars

#KROME_coolingZ_rates

      coolingZ_rates(:) = k(:)

      !check rates > 1
      if(maxval(k)>1d0) then
         print *,"ERROR: found rate >1d0 in coolingZ_rates!"
         print *," Tgas =",Tgas
         do i=1,nZrate
            if(k(i)>1d0) print *,i,k(i)
         end do
         stop
      end if

      !check rates <0
      if(minval(k)<0d0) then
         print *,"ERROR: found rate <0d0 in coolingZ_rates!"
         print *," Tgas =",Tgas
         do i=1,nZrate
            if(k(i)<0d0) print *,i,k(i)
         end do
         stop
      end if

    end function coolingZ_rates

    !**********************
    function coolingZ_rate_tabs(inTgas)
      use krome_commons
      implicit none
      real*8::inTgas,Tgas,coolingZ_rate_tabs(nZrate),k(nZrate)
      integer::idx,j
      Tgas = inTgas

      idx = (log10(Tgas)-coolTab_logTlow) * inv_coolTab_idx + 1

      idx = max(idx,1)
      idx = min(idx,coolTab_n-1)

      do j=1,nZrate
         k(j) = (Tgas-coolTab_T(idx)) * inv_coolTab_T(idx) * &
              (coolTab(j,idx+1)-coolTab(j,idx)) + coolTab(j,idx)
         k(j) = max(k(j), 0d0)
      end do

      coolingZ_rate_tabs(:) = k(:)

    end function coolingZ_rate_tabs

    !**********************
    subroutine coolingZ_init_tabs()
      use krome_commons
      implicit none
      integer::j,jmax,idx
      real*8::Tgas,Tgasold

      jmax = coolTab_n !size of the cooling tables (number of saples)

      !note: change upper and lower limit for rate tables here
      coolTab_logTlow = log10(2d0)
      coolTab_logTup = log10(1d8)

      !pre compute this value since used jmax times
      inv_coolTab_idx = (jmax-1) / (coolTab_logTup-coolTab_logTlow)

      !loop over the jmax interpolation points
      do j=1,jmax
         !compute Tgas for the given point
         Tgas = 1d1**((j-1)*(coolTab_logTup-coolTab_logTlow) &
              /(jmax-1) + coolTab_logTlow)
         !produce cooling rates for the given Tgas
         coolTab(:,j) = coolingZ_rates(Tgas)
         !store Tgas into the array
         coolTab_T(j) = Tgas
         !save 1/dT since it is known
         if(j>1) inv_coolTab_T(j-1) = 1d0 / (Tgas-Tgasold)
         Tgasold = Tgas
      end do

    end subroutine coolingZ_init_tabs

    !*******************************
    !this subroutine solves a non linear system
    ! with the equations stored in fcn function
    ! and a dummy jacobian jcn
    subroutine nleq_wrap(x)
      use krome_user_commons
      integer,parameter::nmax=100 !problem size
      integer,parameter::liwk=nmax+50 !size integer workspace
      integer,parameter::lrwk=(nmax+13)*nmax+60 !real workspace
      integer,parameter::luprt=6 !logical unit verbose output
      integer::neq,iopt(50),ierr,niw,nrw,iwk(liwk),ptype,i
      real*8::x(:),xscal(nmax),rtol,rwk(lrwk),idamp,mdamp,xi(size(x)),minx
      real*8::store_invdvdz
      neq = size(x)
      niw = neq+50
      nrw = (neq+13)*neq+60

      ptype = 2 !initial problem type, 2=mildly non-linear
      rtol = 1d-5 !realtive tolerance
      xi(:) = x(:) !store initial guess
      idamp = 1d-4 !initial damp (when ptype>=4, else default)
      mdamp = 1d-8 !minimum damp (when ptype>=4, else default)
      ierr = 0

      !iterate until ierr==0 and non-negative solutions
      do
         if(ptype>50) then
            print *,"ERROR in nleq1: can't find a solution after attempt",ptype
            stop
         end if

         x(:) = xi(:) !restore initial guess

         !if damping error or negative solutions
         ! prepares initial guess with the thin case
         if(ptype>7.and.(ierr==3.or.ierr==0)) then
            rtol = 1d-5
            iwk(:) = 0
            iopt(:) = 0
            rwk(:) = 0d0
            xscal(:) = 0d0
            store_invdvdz = krome_invdvdz !store global variable
            krome_invdvdz = 0d0 !this sets beta to 1
#IFKROME_use_NLEQ
            call nleq1(neq,fcn,jcn,x(:),xscal(:),rtol,iopt,ierr,&
                 liwk,iwk(:),lrwk,rwk(:))
#ENDIFKROME_use_NLEQ
            if(ierr.ne.0) then
               print *,"ERROR in nleq for thin approx",ierr
               stop
            end if
            krome_invdvdz = store_invdvdz !restore global variable
         end if
         xscal(:) = 0d0 !scaling factor
         rtol = 1d-5 !relative tolerance
         iwk(:) = 0 !default iwk
         iwk(31) = int(1e8) !max iterations
         iopt(:) = 0 !default iopt
         iopt(31) = min(ptype,4) !problem type
         rwk(:) = 0d0 !default rwk
         !reduce damps if damping error
         if(ptype>4.and.ierr==3) then
            idamp = idamp * 1d-1 !reduce idamp
            mdamp = mdamp * 1d-1 !reduce mdamp
         end if
         !if problem is extremely nonlinear use custom damps
         if(ptype>4) then
            rwk(21) = idamp !copy idamp to solver
            rwk(22) = mdamp !copy mdamp to solver
         end if

#IFKROME_use_NLEQ
         call nleq1(neq,fcn,jcn,x(:),xscal(:),rtol,iopt,ierr,&
              liwk,iwk(:),lrwk,rwk(:))
#ENDIFKROME_use_NLEQ

         !check for errors
         if(ierr.ne.0) then
            !print *,"error",ierr
            !problem with damping factor and/or problem type
            if(ierr==3) then
               ptype = ptype + 1 !change the problem type (non-linearity)
            elseif(ierr==5) then
               xi(:) = x(:)
            else
               !other type of error hence stop
               print *,"ERROR in nleq1, ierr:",ierr
               print *,"solutions found so far:"
               do i=1,size(x)
                  print *,i,x(i)
               end do
               stop
            end if
         else
            !if succesful search for negative results
            minx = minval(x) !minimum value
            !if minimum value is positive OK
            if(minx.ge.0d0) then
               exit
            else
               !if negative values are small set to zero
               if(abs(minx)/maxval(x)<rtol) then
                  do i=1,neq
                     x(i) = max(x(i),0d0)
                  end do
                  exit
               else
                  !if large negative values increase non-linearity
                  ptype = ptype + 1
               end if
            end if
         end if
      end do
    end subroutine nleq_wrap

    !***************************
    subroutine fcn(n,x,f,ierr)
      implicit none
      integer::n,ierr
      real*8::x(n),f(n)

#KROME_fcn_cases

    end subroutine fcn

    !**********************************
    !dummy jacobian for non linear equation solver
    subroutine jcn()

    end subroutine jcn

#KROME_coolingZ_functions

#ENDIFKROME

    !***********************
    subroutine mylin2(a,b)
      !solve Ax=B analytically for a 2-levels system
      implicit none
      integer,parameter::n=2
      real*8::a(n,n),b(n),c(n),iab

      !uncomment this: safer but slower function
      !if(a(2,2)==a(2,1)) then
      !   print *,"ERROR: a22=a21 in mylin2"
      !   stop
      !end if
      iab = b(1)/(a(2,2)-a(2,1))
      c(1) = a(2,2) * iab
      c(2) = -a(2,1) * iab
      b(:) = c(:)

    end subroutine mylin2


    !************************
    subroutine mylin3(a,b)
      !solve Ax=B analytically for a 3-levels system
      implicit none
      integer,parameter::n=3
      real*8::iab,a(n,n),b(n),c(n)

      !uncomment this: safer but slower function
      !if(a(2,2)==a(2,3)) then
      !   print *,"ERROR: a22=a23 in mylin3"
      !   stop
      !end if

      !uncomment this: safer but slower
      !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
      !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
      !   print *,"ERROR: division by zero in mylin3"
      !   stop
      !end if

      iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
           + a(2,3)*(a(3,2)-a(3,1)))
      c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
      c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
      c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
      b(:) = c(:)

    end subroutine mylin3

    !************************************
    subroutine plot_cool(n)
      !routine to plot cooling at runtime
      real*8::n(:),Tgas,Tmin,Tmax
      real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
      real*8::cool_dH,cool_Z
      integer::i,imax
      imax = 1000
      Tmin = log10(1d1)
      Tmax = log10(1d8)
      print *,"plotting cooling..."
      open(33,file="KROME_cooling_plot.dat",status="replace")
      do i=1,imax
         Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
         cool_H2 = 0.d0
         cool_H2GP = 0.d0
         cool_HD = 0.d0
         cool_atomic = 0.d0
         cool_Z = 0.d0
         cool_dH = 0.d0
#IFKROME_useCoolingH2
         cool_H2 = cooling_H2(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingH2GP
         cool_H2GP = cooling_H2GP(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingAtomic
         cool_atomic = cooling_atomic(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingHD
         cool_HD = cooling_HD(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingZ
         cool_Z = cooling_Z(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingdH
         cool_dH = cooling_dH(n(:),Tgas)
#ENDIFKROME
         cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
         cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
         write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
              cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
      end do
      close(33)
      print *,"done!"

    end subroutine plot_cool

    !***********************************
    !routine to dump cooling in unit nfile
    subroutine dump_cool(n,Tgas,nfile)
      use krome_commons
      implicit none
      real*8::Tgas,n(:),cools(ncools)
      integer::nfile

      cools(:) = get_cooling_array(n(:),Tgas)
      write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

    end subroutine dump_cool

  end module KROME_cooling
