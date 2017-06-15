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
  subroutine init_photoBins(Tgas)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_getphys
    implicit none
    integer::i,j
    real*8::Tgas,imass(nspec),kt2
    real*8::energy_eV,kk,energyL,energyR,dshift(nmols)

    !rise error if photobins are not defined
    if(photoBinEmid(nPhotoBins)==0d0) then
       print *,"ERROR: when using photo bins you must define"
       print *," the energy interval in bins!"
       stop
    end if

    !get inverse of mass
    imass(:) = get_imass()

    !precompute adimensional line broadening
#KROME_broadening_shift_precalc

#KROME_load_xsecs_from_file

    !tabulate the xsecs into a bin-based array
    do j=1,nPhotoBins
       energyL = photoBinEleft(j)
       energyR = photoBinEright(j)
       energy_eV = photoBinEmid(j) !energy of the bin in eV
#KROME_photobin_xsecs
    end do

    !save interpolated xsecs to file
#KROME_save_xsecs_to_file

    !energy tresholds (eV)
#KROME_photobin_Eth

    !interpolate dust qabs
#KROME_interpolate_dust_qabs

    !map with X->B/C transition to bin corrspondence
#KROME_init_H2kpd_transition_map

  end subroutine init_photoBins

  !**********************
  !save xsecs with index idx to file
  subroutine save_xsec(fname,idx)
    use krome_commons
    implicit none
    character(len=*)::fname
    integer::idx,j
    real*8::energyLeft,energyRight

    open(22,file=trim(fname),status="replace")
    do j=1,nPhotoBins
       energyLeft = photoBinELeft(j) !left bin energy, eV
       energyRight = photoBinERight(j) !right bin energy, eV
       write(22,*) energyLeft, energyRight, photoBinJTab(idx,j)
    end do
    close(22)

  end subroutine save_xsec

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
    use krome_getphys
    implicit none
    integer::i,j
    real*8::dE,kk,Jval,E,Eth,n(:),ncol(nmols),tau

#IFKROME_usePhotoOpacity
    !get column density from number density
    do i=1,nmols
       ncol(i) = num2col(n(i),n(:))
    end do
#ENDIFKROME_usePhotoOpacity

    !init rates and heating
    photoBinRates(:) = 0d0 !1/s/Hz
    photoBinHeats(:) = 0d0 !eV/s/Hz
    GHabing_thin = 0d0 !habing flux
    !loop on energy bins
    do j=1,nPhotoBins
       dE = photoBinEdelta(j) !energy interval, eV
       E = photoBinEmid(j) !energy of the bin in eV
       Jval = photoBinJ(j) !radiation intensity eV/s/cm2/sr/Hz
       if(E>=6d0.and.E<=13.6)then
          GHabing_thin = GHabing_thin + Jval * dE
       endif
       tau = 0d0
#KROME_photobin_opacity
       !loop on reactions
       do i=1,nPhotoRea
          Eth = photoBinEth(i) !reaction energy treshold, eV
          if(E>Eth) then
             !approx bin integral
             kk = photoBinJTab(i,j)*Jval/E*dE
             photoBinRates(i) = photoBinRates(i) + kk
#IFKROME_photobin_heat
             photoBinHeats(i) = photoBinHeats(i) + kk*(E-Eth)
#ENDIFKROME_photobin_heat
          end if
       end do
    end do

    !Final Habing flux
    GHabing_thin = GHabing_thin * 4d0 * pi / (1.6d-3) * iplanck_eV * eV_to_erg

    !converts to 1/s
    photoBinRates(:) = 4d0*pi*photoBinRates(:) * iplanck_eV

#IFKROME_photobin_heat
    !converts to erg/s
    photoBinHeats(:) = 4d0*pi*photoBinHeats(:) * iplanck_eV * eV_to_erg
#ENDIFKROME_photobin_heat

  end subroutine calc_photoBins_thick

  !********************
  !Verner+96 cross section fit (cm2)
  function sigma_v96(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
    implicit none
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

  !************************
  !load the xsecs from file and get limits
  subroutine load_xsec(fname,xsec_val,xsec_Emin,xsec_n,xsec_idE)
    implicit none
    real*8,allocatable::xsec_val(:)
    real*8::xsec_Emin,xsec_dE,xsec_val_tmp(int(1e6)),rout(2)
    real*8::xsec_E_tmp(size(xsec_val_tmp)),xsec_idE,diff
    integer::xsec_n,ios
    character(*)::fname

    !if file already loaded skip subroutine
    if(allocated(xsec_val)) return

    xsec_n = 0 !number of lines found
    !open file
    open(33,file=fname,status="old",iostat=ios)
    !check if file exists
    if(ios.ne.0) then
       print *,"ERROR: problems loading "//fname
       stop
    end if

    !read file line-by-line
    do
       read(33,*,iostat=ios) rout(:) !read line
       if(ios<0) exit !eof
       if(ios/=0) cycle !skip blanks
       xsec_n = xsec_n + 1 !increase line number
       xsec_val_tmp(xsec_n) = rout(2) !read xsec value cm2
       xsec_E_tmp(xsec_n) = rout(1) !read energy value eV
       !compute the dE for the first interval
       if(xsec_n==2) xsec_dE = xsec_E_tmp(2)-xsec_E_tmp(1)
       !check if all the intervals have the same spacing
       if(xsec_n>2) then
          diff = xsec_E_tmp(xsec_n)-xsec_E_tmp(xsec_n-1)
          if(abs(diff/xsec_dE-1d0)>1d-6) then
             print *,"ERROR: spacing problem in file "//fname
             print *," energy points should be equally spaced!"
             print *,"Point number: ",xsec_n
             print *,"Found ",diff
             print *,"Should be",xsec_dE
             stop
          end if
       end if
    end do
    close(33)

    !store the minimum energy
    xsec_Emin = xsec_E_tmp(1)
    !allocate the array with the values
    allocate(xsec_val(xsec_n))
    !copy the values from the temp array to the allocated one
    xsec_val(:) = xsec_val_tmp(1:xsec_n)
    !store the inverse of the delta energy
    xsec_idE = 1d0 / xsec_dE

  end subroutine load_xsec

  !**********************
  !return averaged xsec in the energy range [xL,xR]
  ! units: eV, cm2; broadening shift is adimensional
  function xsec_interp(xL,xR,xsec_val,xsec_Emin,xsec_idE,dshift) result(xsecA)
    use krome_user_commons
    implicit none
    real*8::xsecA,dE,dshift,dE_shift,eL,eR,dxi
    real*8::energy,xsec_val(:),xsec_Emin,xsec_idE,xL,xR
    integer::idx

    !xsec energy step (regular grid)
    dE = 1d0/xsec_idE
    !store inverse of bin size
    dxi = 1d0/(xR-xL)
    xsecA = 0d0 !init integrated xsec
    !loop on xsec vals
    do idx=1,size(xsec_val)
       eL = (idx-1)*dE+xsec_Emin !left interval
       eR = eL + dE !right interval
       energy = (eL+eR)/2d0 !mid point

       !compute line broadening
       eL = eL - 0.5d0*dshift*energy
       eR = eR + 0.5d0*dshift*energy

       !if xsec energy in the interval compute area
       if(xR<eL.and.xL<eL) then
          xsecA = xsecA + 0d0
       elseif(xR>eL.and.xL>eL) then
          xsecA = xsecA + 0d0
       else
          !renormalize xsec area considering partial overlap
          xsecA = xsecA +xsec_val(idx) * (min(eR,xR)-max(eL,xL)) * dxi #KROME_xsecKernelFunction
       end if
    end do

  end function xsec_interp

  !**********************
  !linear interpolation for the photo xsec
  function xsec_interp_mid(energy,xsec_val,xsec_Emin,xsec_n,xsec_idE)
    implicit none
    real*8::xsec_interp_mid,E0
    real*8::energy,xsec_val(:),xsec_Emin,xsec_idE
    integer::xsec_n,idx

    xsec_interp_mid = 0d0
    !retrive index
    idx = (energy-xsec_Emin) * xsec_idE + 1

    !lower bound
    E0 = xsec_Emin + (idx-1)/xsec_idE

    !out of the limits is zero
    if(idx<1.or.idx>xsec_n-1) return

    !linear interpolation
    xsec_interp_mid = (energy-E0) * xsec_idE &
         * (xsec_val(idx+1)-xsec_val(idx)) + xsec_val(idx)

    !avoid negative xsec values when outside the limits
    xsec_interp_mid = max(xsec_interp_mid,0d0)

  end function xsec_interp_mid

  !************************
  !load photodissociation data from default file
  subroutine kpd_H2_loadData()
    use krome_commons
    implicit none
    integer::unit,ios,ii,jj
    real*8::xE,dE,pre
    character(len=20)::fname

    !open file to read
    fname = "H2pdB.dat"
    open(newunit=unit,file=trim(fname),status="old",iostat=ios)
    !check for errors
    if(ios/=0) then
       print *,"ERROR: problem loading file "//trim(fname)
       stop
    end if

    !init data default
    H2pdData_EX(:) = 0d0
    H2pdData_dE(:,:) = 0d0
    H2pdData_pre(:,:) = 0d0

    !loop on file to read
    do
       read(unit,*,iostat=ios) ii,jj,xE,dE,pre
       !skip comments
       if(ios==59.or.ios==5010) cycle
       !exit when eof
       if(ios/=0) exit
       !store data
       H2pdData_EX(ii+1) = xE !ground level energy, eV
       H2pdData_dE(ii+1,jj+1) = dE !Ej-Ei energy, eV
       H2pdData_pre(ii+1,jj+1) = pre !precomp (see file header)
    end do

    !check if enough data have been loaded (file size is expected)
    if((ii+1/=H2pdData_nvibX).or.(jj+1/=H2pdData_nvibB)) then
       !print error message
       print *,"ERROR: missing data when loading "//fname
       print *,"found:",ii+1,jj+1
       print *,"expected:",H2pdData_nvibX,H2pdData_nvibB
       stop
    end if

    close(unit)

  end subroutine kpd_H2_loadData

  !************************
  subroutine kpd_bin_map()
    use krome_commons
    implicit none
    integer::i,j,k
    logical::found

    !loop on excited states (B)
    do i=1,H2pdData_nvibB
       !loop on ground states (X)
       do j=1,H2pdData_nvibX
          !if prefactor is zero no need to check map
          ! default is set to 1 (be aware of it!)
          if(H2pdData_pre(j,i)==0d0) then
             H2pdData_binMap(j,i) = 1
             cycle
          end if

          found = .false.
          !loop on bins
          do k=1,nPhotoBins
             !find energy bin corresponding on the given dE
             if((photoBinEleft(k).le.H2pdData_dE(j,i)) &
                  .and. (photoBinEright(k).ge.H2pdData_dE(j,i))) then
                H2pdData_binMap(j,i) = k
                found = .true.
             end if
          end do
          !error if outside bounds
          if(.not.found) then
             print *,"ERROR: problem when creating H2"
             print *," photodissociation map!"
             print *," min/max (eV):", minval(photoBinEleft), &
                  maxval(photoBinEright)
             print *," transition:",j,i
             print *," corresponding energy (eV):",H2pdData_dE(j,i)
             print *," transitions min/max (eV):", &
                  minval(H2pdData_dE, mask=((H2pdData_dE>0d0) .and. &
                  (H2pdData_pre>0d0))), &
                  maxval(H2pdData_dE, mask=(H2pdData_pre>0d0))
             stop
          end if
       end do
    end do

  end subroutine kpd_bin_map

  !************************
  !compute vibrational partition function at given Tgas
  ! for all the loaded energies (for H2 Solomon)
  function partitionH2_vib(Tgas) result(z)
    use krome_constants
    use krome_commons
    implicit none
    real*8::Tgas,z(H2pdData_nvibX),b
    integer::j

    !prepare partition function from ground (X) levels energies
    b = iboltzmann_eV/Tgas
    z(:) = exp(-H2pdData_EX(:)*b)

    !normalize
    z(:) = z(:)/sum(z)

  end function partitionH2_vib

  !************************
  !compute H2 photodissociation rate (Solomon)
  ! state to state, using preloded data, 1/s
  function kpd_H2(Tgas) result(kpd)
    use krome_commons
    implicit none
    integer::i,j
    real*8::Tgas,kpd,dE,z(H2pdData_nvibX)

    !get partition for ground state X
    z(:) = partitionH2_vib(Tgas)

    !compute the rate, using preloaded data
    kpd = 0d0
    !loop on excited states (B)
    do i=1,H2pdData_nvibB
       !compute rate for ith state
       kpd = kpd + sum(H2pdData_pre(:,i) &
            * photoBinJ(H2pdData_binMap(:,i)) * z(:))
    end do

  end function kpd_H2

  !************************
  !photodissociation H2 xsec from atomic data (for opacity)
  function kpd_H2_xsec(Tgas) result(xsec)
    use krome_constants
    use krome_commons
    implicit none
    real*8::xsec(nPhotoBins),z(H2pdData_nvibX)
    real*8::Tgas
    integer::i

    !get partition for ground state X
    z(:) = partitionH2_vib(Tgas)

    xsec(:) = 0d0
    !loop on excited states (B)
    do i=1,H2pdData_nvibB
       xsec(H2pdData_binMap(:,i)) = &
            xsec(H2pdData_binMap(:,i)) &
            + H2pdData_pre(:,i)*z(:)
    end do

    !cm2
    xsec(:) = xsec(:)*planck_eV

  end function kpd_H2_xsec

  !************************
  !H2 direct photodissociation in the Lyman-Werner bands
  ! cross-section in cm^2 fit by Abel et al. 1997 of
  ! data by Allison&Dalgarno 1969
  function H2_sigmaLW(energy_eV)
    use krome_commons
    implicit none
    real*8::H2_sigmaLW,energy_eV
    real*8::sL0,sW0,sL1,sW1,fact

    !initialization
    sL0 = 0d0
    sL1 = 0d0
    sW0 = 0d0
    sW1 = 0d0

    if(energy_eV>14.675.and.energy_eV<16.820)then
       sL0 = 1d-18*1d1**(15.1289-1.05139*energy_eV)
    elseif(energy_eV>16.820.and.energy_eV<17.6d0)then
       sL0 = 1d-18*1d1**(-31.41d0+1.8042d-2*energy_eV**3-4.2339d-5*energy_eV**5)
    endif

    if(energy_eV>14.675d0.and.energy_eV<17.7d0)then
       sW0 = 1d-18*1d1**(13.5311d0-0.9182618*energy_eV)
    endif

    if(energy_eV>14.159d0.and.energy_eV<15.302d0)then
       sL1 = 1d-18*1d1**(12.0218406d0-0.819429*energy_eV)
    elseif(energy_eV>15.302d0.and.energy_eV<17.2d0)then
       sL1 = 1d-18*1d1**(16.04644d0-1.082438*energy_eV)
    endif

    if(energy_eV>14.159d0.and.energy_eV<17.2d0)then
       sW1 = 1d-18*1d1**(12.87367-0.85088597*energy_eV)
    endif

    fact = 1d0/(phys_orthoParaRatio+1d0)

    H2_sigmaLW = fact*(sL0+sW0)+(1d0-fact)*(sL1+sW1)

  end function H2_sigmaLW
#ENDIFKROME

end module krome_photo
