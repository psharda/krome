module krome_user
  implicit none

#KROME_header

#KROME_species

#KROME_cool_index

#KROME_heat_index

#KROME_common_alias

#KROME_constant_list

contains

#KROME_user_commons_functions

  !***********************
  !print the current values of the phys variables
  subroutine krome_print_phys_variables()
    use krome_commons
    implicit none

#KROME_print_phys_variables

  end subroutine krome_print_phys_variables

#KROME_set_get_phys_functions


#IFKROME_useXrays
  !****************************
  !set the value of J21xrays for tabulated
  ! heating and rate
  subroutine krome_set_J21xray(xarg)
    use krome_commons
    implicit none
    real*8::xarg

    J21xray = xarg

  end subroutine krome_set_J21xray
#ENDIFKROME

#KROME_cooling_functions

#IFKROME_use_coolingZ
  !***************************
  !dump the population of the Z cooling levels
  ! in the nfile file unit, using Tgas as 
  ! independent variable. alias of
  ! dump_cooling_pop subroutine
  subroutine krome_popcool_dump(inTgas,nfile)
    use krome_cooling
    implicit none
    real*8::Tgas,inTgas
    integer::nfile

    Tgas = inTgas
    call dump_cooling_pop(Tgas,nfile)

  end subroutine krome_popcool_dump
#ENDIFKROME

#IFKROME_useDust

  !*************************
  !this subroutine sets the dust distribution in the range
  ! alow_arg, aup_arg, using power law with exponent phi_arg.
  ! All these arguments are optional.
  subroutine krome_set_dust_distribution(alow_arg,aup_arg,phi_arg)
    use krome_dust
    real*8,optional::alow_arg,aup_arg,phi_arg
    real*8::alow,aup,phi

    !default values
    alow = 5d-7 !lower size (cm)
    aup = 2.5d-5 !upper size (cm)
    phi = -3.5d0 !MNR distribution exponent (with its sign)

    if(present(alow_arg)) alow = alow_arg
    if(present(aup_arg)) aup = aup_arg
    if(present(phi_arg)) phi = phi_arg

    call set_dust_distribution(alow,aup,phi)

  end subroutine krome_set_dust_distribution

  !*****************************
  !this function returns an array of size krome_ndust
  ! that contains the amount of dust per bin in 1/cm3.
  function krome_get_dust_distribution()
    use krome_commons
    implicit none
    real*8::krome_get_dust_distribution(ndust)

    krome_get_dust_distribution(:) = xdust(:)

  end function krome_get_dust_distribution

  !******************************
  !this function returns an array of size krome_ndust
  ! that contains the size of the dust bins in cm
  function krome_get_dust_size()
    use krome_commons
    implicit none
    real*8::krome_get_dust_size(ndust)

    krome_get_dust_size(:) = krome_dust_asize(:)

  end function krome_get_dust_size

  !************************
  !this function sets the default temperature
  ! for all the dust bins.
  subroutine krome_set_defaultTdust(arg)
    use krome_commons
    implicit none
    real*8::arg

    krome_dust_T(:) = arg

  end subroutine krome_set_defaultTdust

  !****************************
  subroutine krome_scale_dust_distribution(xscale)
    use krome_commons
    implicit none
    real*8::xscale

    xdust(:) = xdust(:) * xscale

  end subroutine krome_scale_dust_distribution

  !***************************
  subroutine krome_scale_dust_gas_ratio(dust_to_gas_ratio,x)
    use krome_commons
    use krome_subs
    implicit none
    real*8::dust_to_gas_ratio,x(:),rho_gas,m(nspec)
    integer::nd

    nd = ndust / ndustTypes
    m(:) = get_mass()
    rho_gas = sum(x(:)*m(1:nmols))
    xdust(:) = dust_to_gas_ratio * rho_gas / krome_grain_rho &
         / krome_dust_asize3(:) / nd
    
  end subroutine krome_scale_dust_gas_ratio

#ENDIFKROME

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on()
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off()
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

#IFKROME_usePhotoBins
  !************************
  subroutine krome_calc_photobins()
    use krome_photo
    call calc_photobins()
  end subroutine krome_calc_photobins

  !****************************
  ! set the energy per photo bin
  ! eV/s/cm2/sr/Hz
  subroutine krome_set_photoBinJ(phbin)
    use krome_commons
    use krome_photo
    implicit none
    real*8::phbin(:)
    photoBinJ(:) = phbin(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBinJ

  !*************************
  ! set the energy (frequency) of the photobin
  ! as left-right limits in eV
  subroutine krome_set_photobinE_lr(phbinleft,phbinright)
    use krome_commons
    use krome_photo
    implicit none
    real*8::phbinleft(:),phbinright(:)
    photoBinEleft(:) = phbinleft(:)
    photoBinEright(:) = phbinright(:)
    photoBinEmid(:) = 0.5d0*(phbinleft(:)+phbinright(:))
    photoBinEdelta(:) = phbinright(:)-phbinleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_lr

  !********************************
  ! set the energy (eV) of the photobin
  ! linearly from lowest and highest energy value
  subroutine krome_set_photobinE_lin(lower,upper)
    use krome_commons
    use krome_photo
    implicit none
    real*8::lower,upper,dE
    integer::i
    dE = abs(upper-lower)/nPhotoBins
    do i=1,nPhotoBins
       photoBinEleft(i) = dE*(i-1) + lower
       photoBinEright(i) = dE*i + lower
       photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
    end do
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_lin

  !********************************
  ! set the energy (eV) of the photobin
  ! logarithmic from lowest to highest energy value
  subroutine krome_set_photobinE_log(lower,upper)
    use krome_commons
    use krome_photo
    implicit none
    real*8::lower,upper,dE,logup,loglow
    integer::i
    if(lower.ge.upper) then
       print *,"ERROR: in  krome_set_photobinE_log lower >= upper limit!"
       stop
    end if
    loglow = log10(lower)
    logup = log10(upper)
    dE = 1d1**(abs(logup-loglow)/nPhotoBins)
    do i=1,nPhotoBins
       photoBinEleft(i) = 1d1**((i-1)*(logup-loglow)/nPhotoBins + loglow)
       photoBinEright(i) = 1d1**(i*(logup-loglow)/nPhotoBins + loglow)
       photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
    end do
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)

    !initialize xsecs table
    call init_photoBins()

  end subroutine krome_set_photobinE_log

  !*********************************
  function krome_get_photoBinJ()
    !returns an array containing the flux for each photo bin
    use krome_commons
    real*8::krome_get_photoBinJ(nPhotoBins)
    krome_get_photoBinJ(:) = photoBinJ(:)
  end function krome_get_photoBinJ

  !*********************************
  function krome_get_photoBinE_left()
    !returns an array with the left energy limits (eV)
    use krome_commons
    real*8::krome_get_photoBinE_left(nPhotoBins)
    krome_get_photoBinE_left(:) = photoBinEleft(:)
  end function krome_get_photoBinE_left

  !*********************************
  function krome_get_photoBinE_right()
    !returns an array with the right energy limits (eV)
    use krome_commons
    real*8::krome_get_photoBinE_right(nPhotoBins)
    krome_get_photoBinE_right(:) = photoBinEright(:)
  end function krome_get_photoBinE_right

  !*********************************
  function krome_get_photoBinE_mid()
    !returns an array with the middle energy limits (eV)
    use krome_commons
    real*8::krome_get_photoBinE_mid(nPhotoBins)
    krome_get_photoBinE_mid(:) = photoBinEmid(:)
  end function krome_get_photoBinE_mid

  !*********************************
  function krome_get_photoBinE_delta()
    !returns an array with the middle energy limits (eV)
    use krome_commons
    real*8::krome_get_photoBinE_delta(nPhotoBins)
    krome_get_photoBinE_delta(:) = photoBinEdelta(:)
  end function krome_get_photoBinE_delta

  !*********************************
  function krome_get_photoBinE_idelta()
    !returns an array with the middle energy limits (eV)
    use krome_commons
    real*8::krome_get_photoBinE_idelta(nPhotoBins)
    krome_get_photoBinE_idelta(:) = photoBinEidelta(:)
  end function krome_get_photoBinE_idelta

  !*********************************
  function krome_get_photoBin_rates()
    !returns an array with the integrated photo rates (1/s)
    use krome_commons
    real*8::krome_get_photoBin_rates(nPhotoRea)
    krome_get_photoBin_rates(:) = photoBinRates(:)
  end function krome_get_photoBin_rates

  !*********************************
  !returns an array with the integrated photo heatings (erg/s)
  function krome_get_photoBin_heats()
    use krome_commons
    real*8::krome_get_photoBin_heats(nPhotoRea)
    krome_get_photoBin_heats(:) = photoBinHeats(:)
  end function krome_get_photoBin_heats

  !****************************
  !multiply all the bins by a factor xscale
  subroutine krome_photoBin_scale(xscale)
    use krome_commons
    use krome_photo
    implicit none
    real*8::xscale

    photoBinJ(:) = photoBinJ(:) * xscale

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale

  !****************************
  !multiply all the bins by an array xscale(:)
  ! of size krome_nPhotoBins
  subroutine krome_photoBin_scale_array(xscale)
    use krome_commons
    use krome_photo
    implicit none
    real*8::xscale(:)

    photoBinJ(:) = photoBinJ(:) * xscale(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale_array

  !**********************************
  subroutine krome_set_photoBin_BBlog(lower,upper,Tbb)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    implicit none
    real*8::lower,upper,Tbb,x,xmax,xexp,Jlim
    integer::i

    !limit for the black body intensity to check limits
    Jlim = 1d-3

    call krome_set_photoBinE_log(lower,upper)

    !eV/cm2/s/Hz/sr
    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       photoBinJ(i) = planckBB(x,Tbb)
    end do

    !uncomment this below for additional control
!!$    !find the maximum using Wien's displacement law
!!$    xmax = Tbb/2.8977721d-1 * clight * planck_eV !eV
!!$
!!$    if(xmax<lower) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is below the lowest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"lowest (eV)",lower
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(xmax>upper) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is above the highest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"highest (eV)",upper
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(1)>Jlim) then
!!$       print *,"WARNING: lower bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(1)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(nPhotoBins)>Jlim) then
!!$       print *,"WARNING: upper bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(nPhotoBins)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_BBlog

  !*************************************
  !set the BB spectrum and the limits using bisection
  subroutine krome_set_photoBin_BBlog_auto(Tbb)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::Tbb,xlow,xup,eps,xmax,J0,J1,x0,x1,xm,Jm
    eps = 1d-6

    !Rayleigh–Jeans approximation for the minimum energy
    xlow = planck_eV*clight*sqrt(.5d0/Tbb/boltzmann_eV*eps)

    !find energy of the Wien maximum (eV)
    xmax = Tbb / 2.8977721d-1 * clight * planck_eV

    !bisection to find the maximum
    x0 = xmax
    x1 = 2.9d2*Tbb*boltzmann_eV 
    J0 = planckBB(x0,Tbb) - eps 
    J1 = planckBB(x1,Tbb) - eps
    if(J0<0d0.or.J1>0d0) then
       print *,"ERROR: problems with auto planck bisection!"
       stop
    end if

    do 
       xm = 0.5d0*(x0+x1)
       Jm = planckBB(xm,Tbb) - eps
       if(Jm>0d0) x0 = xm
       if(Jm<0d0) x1 = xm
       if(abs(Jm)<eps*1d-3) exit
    end do
    xup = xm

    !initialize BB radiation using the values found
    call krome_set_photoBin_BBlog(xlow,xup,Tbb)


  end subroutine krome_set_photoBin_BBlog_auto

  !**************************
  subroutine krome_set_photoBin_draineLin(lower,upper)
    use krome_commons
    use krome_photo
    use krome_constants
    real*8::upper,lower,x
    integer::i

    call krome_set_photoBinE_lin(lower,upper)

    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       !eV/cm2/sr/s/Hz
       if(x<13.6d0) then
          photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
               * x *planck_eV
       else
          photoBinJ(i) = 0d0
       end if
    end do

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLin


  !**************************
  subroutine krome_set_photoBin_draineLog(lower,upper)
    use krome_commons
    use krome_photo
    use krome_constants
    real*8::upper,lower,x
    integer::i

    call krome_set_photoBinE_log(lower,upper)

    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       !eV/cm2/sr/s/Hz
       if(x<13.6d0) then
          photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
               * x *planck_eV
       else
          photoBinJ(i) = 0d0
       end if
    end do

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLog

  !**************************
  subroutine krome_set_photoBin_J21lin(lower,upper)
    use krome_commons
    use krome_photo
    real*8::upper,lower

    call krome_set_photoBinE_lin(lower,upper)
    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/s/Hz/sr

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21lin

  !**************************
  subroutine krome_set_photoBin_J21log(lower,upper)
    use krome_commons
    use krome_photo
    real*8::upper,lower

    call krome_set_photoBinE_log(lower,upper)
    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/s/Hz/sr

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21log

  !*****************************
  !get the opacity exp(-tau) correpsonding the to x(:)
  ! chemical composition. The column density
  ! is computed using the expression in the 
  ! num2col(x) function.
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity(x,Tgas)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    real*8::x(:),tau,krome_get_opacity(nPhotoBins),Tgas,n(nspec)
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
       tau = 0d0
       !loop on species
       do i=1,nPhotoRea
          !calc opacity as column_density * cross_section
          idx = photoPartners(i)
          tau = tau + num2col(x(idx),n(:)) * photoBinJTab(i,j)
       end do
       krome_get_opacity(j) = tau !store
    end do

  end function krome_get_opacity

  !*****************************
  !get the opacity exp(-tau) correpsonding to the x(:)
  ! chemical composition. The column density
  ! is computed using the size of the cell (csize)
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity_size(x,Tgas,csize)
    use krome_commons
    use krome_photo
    use krome_subs
    implicit none
    real*8::x(:),tau,krome_get_opacity_size(nPhotoBins),Tgas
    real*8::csize,n(nspec)
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
       tau = 0d0
       !loop on species
       do i=1,nPhotoRea
          !calc opacity as column_density * cross_section
          !where column_density is density*cell_size
          idx = photoPartners(i)
          tau = tau + x(idx) * csize * photoBinJTab(i,j)
       end do
       krome_get_opacity_size(j) = tau !store
    end do

  end function krome_get_opacity_size

  !*******************************
  !dump the Jflux profile to the file
  ! with unit number nfile
  subroutine krome_dump_Jflux(nfile)
    use krome_commons
    implicit none
    integer::i,nfile

    do i=1,nPhotoBins
       write(nfile,*) photoBinEmid(i),photoBinJ(i)
    end do

  end subroutine krome_dump_Jflux

#ENDIFKROME

  !***************************
  !alias for coe in krome_subs
  ! returns the coefficients for a given Tgas
  function krome_get_coef(Tgas)
    use krome_commons
    use krome_subs
    real*8::krome_get_coef(nrea),Tgas,n(nspec)
    n(:) = 0d0
    n(idx_Tgas) = Tgas

    krome_get_coef(:) = coe(n(:))

  end function krome_get_coef

  !****************************
  !get the mean molecular weight from
  ! mass fractions
  function krome_get_mu_x(xin)
    use krome_commons
    implicit none
    real*8::krome_get_mu_x,xin(:)
    real*8::n(nmols)
    n(:) = krome_x2n(xin(:),1d0)
    krome_get_mu_x = krome_get_mu(n(:))
  end function krome_get_mu_x

  !****************************
  !return the adiabatic index from mass fractions
  ! and temperature in K
  function krome_get_gamma_x(xin,inTgas)
    use krome_commons
    implicit none
    real*8::krome_get_gamma_x,rhogas,xin(:)
    real*8::x(nmols),Tgas,inTgas

    Tgas = inTgas
    x(:) = krome_x2n(xin(:),1d0)
    krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

  end function krome_get_gamma_x


  !***************************
  !normalize mass fractions and
  ! set charge to zero
  subroutine krome_consistent_x(x)
    use krome_commons
    use krome_constants
    implicit none
    real*8::x(nmols),isumx,sumx,xerr,imass(nmols),ee

    !1. charge consistency
    imass(:) = krome_get_imass()

#KROME_zero_electrons

    ee = sum(krome_get_charges()*x(:)*imass(:))
    ee = max(ee*e_mass,0d0)
#KROME_electrons_balance

    !2. mass fraction consistency
    sumx = sum(x)

    !NOTE: uncomment here if you want some additional control
    !conservation error threshold: rise an error if above xerr
    !xerr = 1d-2
    !if(abs(sum-1d0)>xerr) then
    !   print *,"ERROR: some problem with conservation!"
    !   print *,"|sum(x)-1|=",abs(sum-1d0)
    !   stop
    !end if

    isumx = 1d0/sumx
    x(:) = x(:) * isumx

  end subroutine krome_consistent_x

  !*********************
  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number 
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    use krome_commons
    implicit none
    real*8::n(nmols),rhogas,krome_n2x(nmols)

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass 
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    use krome_commons
    implicit none
    real*8::x(nmols),rhogas,krome_x2n(nmols)

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt)
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    implicit none
    real*8::x(:),Tgas,dt,dTgas,k(nrea),krome_gamma
    real*8::n(nspec),nH2dust

#IFKROME_use_thermo
    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
         * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas
#ENDIFKROME

  end subroutine krome_thermo

#IFKROME_use_heating
  !*************************
  function krome_get_heating(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::krome_get_heating,x(nmols),n(nspec)
    real*8::Tgas,inTgas,k(nrea),nH2dust
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    k(:) = coe(n(:))
    nH2dust = 0d0
    krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
  end function krome_get_heating

  !*****************************
  function krome_get_heating_array(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::x(:),n(nspec),inTgas,k(nrea)
    real*8::krome_get_heating_array(8),Tgas,nH2dust

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    k(:) = coe(n(:))
    Tgas = inTgas
    nH2dust = 0d0
    krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

  end function krome_get_heating_array
#ENDIFKROME

#IFKROME_use_cooling
  !*************************
  function krome_get_cooling(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::krome_get_cooling,x(nmols),n(nspec)
    real*8::Tgas,inTgas
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    krome_get_cooling = cooling(n,Tgas)
  end function krome_get_cooling

  !*****************************
  function krome_get_cooling_array(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::x(:),n(nspec),inTgas
    real*8::krome_get_cooling_array(ncools),Tgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    Tgas = inTgas
    krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

  end function krome_get_cooling_array

  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n)
    use krome_cooling
    implicit none
    real*8::n(:)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nmols),Tgas,x(nspec)
    integer,optional::nfile_in
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
    if(present(nfile_in)) nfile = nfile_in
    call dump_cool(x(:),Tgas,nfile)

  end subroutine krome_dump_cooling

#ENDIFKROME

  !*************************
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8::x(:),xi(:)
    real*8::n(krome_nspec),ni(krome_nspec),krome_conserve(krome_nmols)

    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas)
    use krome_subs
    use krome_commons
    real*8::krome_get_gamma,x(nmols),n(nspec),Tgas
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies
  !alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    implicit none
    integer::krome_get_zatoms(nmols),zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)

  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from 
  ! number density and mass density
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8::krome_get_mu,x(:),n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    krome_get_mu = get_mu(n(:))
  end function krome_get_mu

  !*****************
  !get an array of double containing the masses in g
  ! of the species
  !alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec), krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse 
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec), krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x)
    use krome_commons
    use krome_subs
    real*8::n(nspec),x(:),krome_get_Hnuclei
    n(:) = 0d0
    n(1:nmols) = x(:)

    krome_get_Hnuclei = get_Hnuclei(n(:))

  end function krome_get_Hnuclei

  !*****************
  !get an array containing the charges of the species
  !alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    implicit none
    real*8::tmp(nspec),krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 containing the names
  ! of alla the species
  !alias for get_names
  function krome_get_names()
    use krome_subs
    use krome_commons
    implicit none
    character*16::krome_get_names(nmols),tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !*****************
  !get the index of the species with name name
  !alias for get_index
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer::krome_get_index
    character*(*)::name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n)
    use krome_commons
    real*8::krome_get_rho,n(:)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho

  !*************************
  subroutine krome_scale_Z(n,Z)
    !scale metallicity the metals contained in n(:) 
    ! to Z according to Arnett(1996)
    use krome_commons
    real*8::n(:),Z,Htot

#KROME_scaleZ

  end subroutine krome_scale_Z

  !***********************
  function krome_get_electrons(n)
    !get the total number of electrons from
    ! the number denisities of all the species
    use krome_commons
    real*8::n(:),ee,krome_get_electrons,x(size(n))
    x(:) = n(:)
#KROME_zero_electrons
    ee = sum(x(:) * krome_get_charges())
    krome_get_electrons = max(0.d0, ee)
  end function krome_get_electrons

  !**********************
  !print the nbest fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest)
    use krome_subs
    use krome_commons
    implicit none
    real*8::x(nmols),xin(nmols),n(nspec),Tgas
    integer::nbest
    x(:) = xin(:)
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)

  end subroutine krome_print_best_flux

  !print the best fluxes greater than a fraction frac
  ! of the maximum flux
  !*********************
  subroutine krome_print_best_flux_frac(xin,Tgas,frac)
    use krome_subs
    use krome_commons
    implicit none
    real*8::x(nmols),xin(nmols),n(nspec),Tgas,frac
    x(:) = xin(:)
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_frac(n,Tgas,frac)

  end subroutine krome_print_best_flux_frac
  
  !**********************
  !print the nbest fluxes for a given species
  subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
    use krome_subs
    use krome_commons
    implicit none
    real*8::x(nmols),xin(nmols),n(nspec),Tgas
    integer::nbest,idx_find
    x(:) = xin(:)
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_spec(n,Tgas,nbest,idx_find)
  end subroutine krome_print_best_flux_spec

  !*******************************
  !get the fluxes of all the reactions in 1/cm3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8::krome_get_flux(nrea),n(nmols),Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*********************
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    implicit none
    real*8::krome_get_qeff(nrea)

    krome_get_qeff(:) = get_qeff()

  end function krome_get_qeff

#IFKROME_useStars

  !**************************
  !alias for stars_coe in krome_star
  function krome_stars_coe(x,rho,Tgas,y_in,zz_in)
    use krome_commons
    use krome_stars
    use krome_subs
    implicit none
    real*8::rho,Tgas,krome_stars_coe(nrea)
    real*8::n(nspec),x(:)
    integer::ny
    real*8,allocatable::y(:)
    integer,allocatable::zz(:)
    real*8,optional::y_in(:)
    integer,optional::zz_in(:)

    !check if extened abundances and zatom array are present
    if(present(y_in)) then
       ny = size(y_in)
       allocate(y(ny), zz(ny))
       y(:) = y_in(:)
       zz(:) = zz_in(:)
    else
       allocate(y(nmols),zz(nmols))
       zz(:) = get_zatoms()
       y(:) = x(:)
    end if

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_stars_coe(:) = stars_coe(n(:),rho,Tgas,y(:),zz(:))
    deallocate(y,zz)
  end function krome_stars_coe

  !********************************
  !alias for stars_energy
  function krome_stars_energy(x,rho,Tgas,k)
    use krome_commons
    use krome_stars
    implicit none
    real*8::x(:),rho,Tgas,krome_stars_energy(nrea)
    real*8::n(nspec),k(nrea)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_stars_energy(:) = stars_energy(n(:),rho,Tgas,k(:))

  end function krome_stars_energy

#ENDIFKROME

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile)
    use krome_commons
    real*8::n(:),Tgas,flux(nrea)
    integer::nfile,i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
       write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
    use krome_commons
    use krome_subs
    implicit none
    integer::funit,i,imax,j
    real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec),inTmin,inTmax

    Tmin = log10(inTmin)
    Tmax = log10(inTmax)

    n(:) = 1d-40
    do i=1,imax
       Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
       n(idx_Tgas) = Tgas
       k(:) = coe(n(:))
       do j=1,nrea
          write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
       end do
       write(funit,*)
    end do


  end subroutine krome_dump_rates

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas)
    use krome_commons
    use krome_subs
    integer::i,charges(nspec)
    real*8::x(:),Tgas,masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
       print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info

end module krome_user
