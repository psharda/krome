module krome_photo
contains

#IFKROME_usePhotoBins

  !*************************
  !get the intensity of the photon flux at
  ! a given energy in eV.
  ! returned value is in eV/cm2/s/Hz
  function get_photoIntensity(energy)
    use krome_commons
    implicit none
    real*8::get_photoIntensity,energy
    integer::i

    !check if requested energy is lower than the lowest limit
    if(energy<photoBinEleft(1)) then
       get_photoIntensity = 0d0 !photoBinJ(1)
       return
    end if
    
    !check if requested energy is greater that the the largest limit
    if(energy>photoBinEright(nPhotoBins)) then
       get_photoIntensity = 0d0 !photoBinJ(nPhotoBins)
       return
    end if

    !look for the interval
    do i=1,nPhotoBins
       if(photoBinEleft(i).le.energy .and. photoBinEright(i).ge.energy) then
          get_photoIntensity = photoBinJ(i)
          return
       end if
    end do

    !error if nothing found
    print *,"ERROR: no interval found in get_photoIntensity"
    print *,"energy:",energy,"eV"
    stop !halt program

  end function get_photoIntensity

  !*********************
  !initialize/tabulate the bin-based xsecs
  subroutine init_photoBins()
    use krome_commons
    implicit none
    integer::i,j
    real*8::energy_eV,kk

    if(photoBinEmid(nPhotoBins)==0d0) then
       print *,"ERROR: when using photo bins you must define"
       print *," the energy interval in bins!"
       stop
    end if

    !tabulate the xsecs into a bin-based array
    do j=1,nPhotoBins
       energy_eV = photoBinEmid(j) !energy of the bin in eV
#KROME_photobin_xsecs
    end do

    !energy tresholds (eV)
#KROME_photobin_Eth

  end subroutine init_photoBins

  !**********************
  !compute integrals to derive phtorates (thin)
  subroutine calc_photoBins()
    use krome_commons
    implicit none
    real*8::n(nspec)
    
    n(:) = 0d0
    call calc_photoBins_thick(n)
    
  end subroutine calc_photoBins

  !**********************
  !compute integrals to derive phtorates (thick)
  subroutine calc_photoBins_thick(n)
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    integer::i,j
    real*8::dE,kk,Jval,E,Eth,n(:),ncol(nmols),tau
    
    !get column density from number density
    do i=1,nmols
       ncol(i) = num2col(n(i))
    end do
    
    !init rates and heating
    photoBinRates(:) = 0d0 !1/s/Hz
    photoBinHeats(:) = 0d0 !eV/s/Hz
    !loop on energy bins
    do j=1,nPhotoBins
       dE = photoBinEdelta(j) !energy interval, eV
       E = photoBinEmid(j) !energy of the bin in eV
       Jval = photoBinJ(j) !radiation intensity eV/s/cm2/sr/Hz
       tau = 0d0
#KROME_photobin_opacity
       !loop on reactions
       do i=1,nPhotoRea
          Eth = photoBinEth(i) !reaction energy treshold, eV
          if(E>Eth) then
             !approx bin integral
             kk = 4d0*pi*photoBinJTab(i,j)*Jval/E*dE * exp(-tau) 
             photoBinRates(i) = photoBinRates(i) + kk
#IFKROME_photobin_heat
             photoBinHeats(i) = photoBinHeats(i) + kk*(E-Eth)
#ENDIFKROME_photobin_heat
          end if
       end do
    end do

    !converts to 1/s
    photoBinRates(:) = 4d0*pi*photoBinRates(:) * iplanck_eV

#IFKROME_photobin_heat
    !converts to erg/s
    photoBinHeats(:) = photoBinHeats(:) * iplanck_eV * eV_to_erg
#ENDIFKROME_photobin_heat

  end subroutine calc_photoBins_thick

#ENDIFKROME  

  !********************
  function sigma_v96(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
    !Verner+96 cross section fit (cm2)
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
    !Heating with Verner+96 cross section fit (cm2*eV)
    use krome_constants
    real*8::heat_v96,energy_eV,sigma_0,Fy,yw,x,y,E0,Eth
    real*8::y0,y1,ya,P
    x = energy_eV/E0 - y0
    y = sqrt(x**2 + y1**2)
    Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
         * (1.d0+sqrt(y/ya))**(-P)
    heat_v96 = 1d-18 * sigma_0 * Fy * (energy_eV - Eth) !cm2*eV
  end function heat_v96

end module krome_photo
