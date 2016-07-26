!This module contains functions and subroutines 
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.
module krome_grfuncs
contains

#KROME_header

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
         * pre_kvgas_sqrt * ims * sqrTgas &
         * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
         + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
       expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
         + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
       dust_get_inv_phi(i) = 0d0
       dephi = (4d0 * nndust(i) * pi * asize2(i))
       if(dephi.le.0d0) cycle
       dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

#IFKROME_useChemisorption
  !***************************
  function dust_get_rateChem_PC(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_PC(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_PC(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_PC(idx+1)-dust_rateChem_PC(idx)) &
            + dust_rateChem_PC(idx)
    end do

  end function dust_get_rateChem_PC

  !***************************
  function dust_get_rateChem_CP(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_CP(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_CP(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_CP(idx+1)-dust_rateChem_CP(idx)) &
            + dust_rateChem_CP(idx)
    end do

  end function dust_get_rateChem_CP

  !***************************
  function dust_get_rateChem_CC(Tdust)
    use krome_commons
    implicit none
    real*8::dust_get_rateChem_CC(ndust), Tdust(ndust)
    integer::i,idx

    do i=1,ndust
       idx = (Tdust(i) - dust_rateChem_xmin) * dust_rateChem_xfact + 1
       dust_get_rateChem_CC(i) = (Tdust(i)-dust_rateChem_x(idx)) * dust_rateChem_invdx &
            * (dust_rateChem_CC(idx+1)-dust_rateChem_CC(idx)) &
            + dust_rateChem_CC(idx)
    end do

  end function dust_get_rateChem_CC
#ENDIFKROME

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
       Td100 = Tdust(i) * 1d-2
       dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
            + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
       dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

#KROME_Ebareice

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

#KROME_Ebareice23

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

#KROME_Ebind_ice

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

#KROME_Ebind_bare

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

#KROME_parent_dust_bin

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

end module krome_grfuncs
