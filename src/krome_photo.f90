module krome_photo
contains

#IFKROME_usePhotoBins

  !*******************************
  !load a frequency-dependent opacity table stored in fname file,
  ! column 1 is energy or wavelenght in un units of unitEnergy
  ! (default eV), column 2 is opacity in cm2/g.
  ! opacity is interpolated over the current photo-binning.
  subroutine load_opacity_table(fname, unitEnergy)
    use krome_commons
    use krome_constants
    implicit none
    integer,parameter::ntmp=int(1e5)
    character(len=*)::fname
    character(len=*),optional::unitEnergy
    character*10::eunit
    integer::ios,icount,iR,iL,i,j,fileUnit
    real*8::wl,opac,fL,fR,kk,dE
    real*8::wls(ntmp),opacs(ntmp)
    real*8,allocatable::energy(:),kappa(:)
    
    !read energy unit optional argument
    eunit = "eV" !default is eV
    if(present(unitEnergy)) then
      eunit = trim(unitEnergy)
    end if
    
    !read form file
    open(newunit=fileUnit,file=trim(fname),status="old",iostat=ios)
    !error if problems reading file
    if(ios/=0) then
      print *,"ERROR: problem while loading "//trim(fname)
      stop
    end if
    icount = 0
    !loop on file lines
    do
      !read wavelength and opacity
      read(fileUnit,*,iostat=ios) wl,opac
      if(ios/=0) exit
      icount = icount + 1
      wls(icount) = wl
      opacs(icount) = opac
    end do
    close(fileUnit)
    
    !allocate arrays
    allocate(energy(icount), kappa(icount))
    !copy temp arrays into allocated arrays, converting units
    if(trim(eunit)=="eV") then
      !eV->eV (default)
      kappa(:) = opacs(1:icount)
      energy(:) = wls(1:icount)
    elseif(trim(eunit)=="micron") then
      !micron->eV
      kappa(:) = opacs(1:icount)
      energy(:) = planck_eV*clight/(wls(1:icount)*1d-4)
    else
      print *,"ERROR: in load opacity table energy unit unknow",trim(eunit)
      stop
    end if

    !reverse array if necessary
    if(energy(2)<energy(1)) then
      energy(:) = energy(size(energy):1:-1)
      kappa(:) = kappa(size(kappa):1:-1)
    end if

    !check if photobins are intialized
    if(maxval(photoBinEleft)==0d0) then
      print *,"ERROR: empty photobins when interpolating dust Qabs"
      print *," from file "//trim(fname)
      print *,"You probably need to define a photobins metric before"
      print *," the call to krome_load_opacity_table"
      stop
    end if

    !check lower limit
    if(photoBinEleft(1)<energy(1)) then
      print *,"ERROR: dust table "//trim(fname)//" energy lower bound (eV)"
      print *,photoBinEleft(1), "<", energy(1)
      stop
    end if

    !check upper limit
    if(photoBinEright(nPhotoBins)>energy(size(energy))) then
      print *,"ERROR: dust table "//trim(fname)//" energy upper bound (eV)"
      print *,photoBinEright(nPhotoBins), ">", energy(size(energy))
      stop
    end if

    !interpolate on current energy distribution
    do j=1,nPhotoBins
      do i=2,size(energy)
        !find left bound position
        if(photoBinEleft(j)>energy(i-1) &
            .and. photoBinEleft(j)<energy(i)) then
        dE = energy(i)-energy(i-1)
        fL = (photoBinEleft(j)-energy(i-1))/dE &
            * (kappa(i)-kappa(i-1)) + kappa(i-1)
        iL = i
      end if

      !find right bound position
      if(photoBinEright(j)>energy(i-1) &
          .and. photoBinEright(j)<energy(i)) then
      dE = energy(i)-energy(i-1)
      fR = (photoBinEright(j)-energy(i-1))/dE &
          * (kappa(i)-kappa(i-1)) + kappa(i-1)
      iR = i
    end if
    end do

    !sum opacity for the given photo bin
    kk = 0d0
    !if there are other opacity points in between left and right limits
    if(iR-iL>0) then
      kk = kk + (energy(iL)-photoBinEleft(j))*(fL+kappa(iL))/2d0
      kk = kk + (photoBinEright(j)-energy(iR-1))*(fR+kappa(iR-1))/2d0
      !sum points in between
      do i=iL,iR-2
        kk = kk + (energy(i+1)-energy(i))*(kappa(i+1)+kappa(i))/2d0
      end do
    elseif(iR==iL) then
      !no opacity points in between
      kk = kk + (fL+fR)*(photoBinEright(j)-photoBinEleft(j))/2d0
    else
      print *,"ERROR: dust opacity interpolation error, iR-iL<0!"
      print *,"iR,iL:",iR,iL
      stop
    end if

    !copy to common and scale to bin size
    dE = photoBinEright(j)-photoBinEleft(j)
    opacityDust(j) = kk/dE

    end do

    !dump interpolated opacity
    open(newunit=fileUnit,file="opacityDust.interp",status="replace")
    do j=1,nPhotoBins
      write(fileUnit,*) photoBinEmid(j),opacityDust(j)
    end do
    close(fileUnit)

    !dump original opacity file (as loaded by krome)
    open(newunit=fileUnit,file="opacityDust.org",status="replace")
    do i=1,size(energy)
      write(fileUnit,*) energy(i),kappa(i)
    end do
    close(fileUnit)

  end subroutine load_opacity_table

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
    dshift(:) = 0d0
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
  !Verner+96 cross section fit (cm2)
  !Average by numerical integration
  function sigma_v96_int(E_low,E_high,E0,sigma_0,ya,P,yw,y0,y1)
    real*8::sigma_v96_int,E_low,E_high,sigma_0,integral,yw,x,y,E0
    real*8::y0,y1,ya,P
    real*8::binWidth,dE,E
    integer::i
    integer,parameter::N=100
    integral = 0d0
    binWidth = E_high-E_low
    dE = binWidth/real(N,kind=8)
    do i=1,N
       E = E_low + (i-0.5)*dE
       integral = integral + sigma_v96(E,E0,sigma_0,ya,P,yw,y0,y1)*dE
    end do
    sigma_v96_int = integral / binWidth !cm2
  end function sigma_v96_int

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
          xsecA = xsecA +xsec_val(idx) * (min(eR,xR)-max(eL,xL)) &
               * dxi #KROME_xsecKernelFunction
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

  ! *****************************
  ! load kabs from file
  subroutine find_Av_load_kabs2(file_name)
    use krome_commons
    use krome_constants
    implicit none
    integer,parameter::imax=10000
    character(len=*),intent(in),optional::file_name
    character(len=200)::fname
    integer::ios, unit, icount, i, j
    real*8::tmp_energy(imax), tmp_data(imax), f1, f2, kavg, ksum
    real*8,allocatable::Jdraine(:)

    ! check if energy bins are set
    if(maxval(photoBinEleft)==0d0) then
       print *, "ERROR: to load kabs for Av G0 finder you"
       print *, " have to initialize some energy bins!"
       stop
    end if

    ! check if optional argument is present
    fname = "kabs_draine_Rv31.dat"
    if(present(file_name)) then
       fname = trim(file_name)
    end if

    ! open file to read
    open(newunit=unit, file=fname, status="old", iostat=ios)
    ! check if file is there
    if(ios/=0) then
       print *, "ERROR: Kabs file not found!"
       print *, trim(fname)
       stop
    end if

    ! loop on file lines
    icount = 1
    do
       read(unit, *, iostat=ios) tmp_energy(icount), &
            tmp_data(icount)
       if(ios/=0) exit
       icount = icount + 1
    end do
    close(unit)

    ! convert microns to eV
    tmp_energy(1:icount-1) = planck_eV * clight &
         / (tmp_energy(1:icount-1) * 1d-4)

    ! get corresponding draine flux
    allocate(Jdraine(icount-1))
    Jdraine(:) = get_draine(tmp_energy(1:icount-1))

    ! loop on photobins to get average kabs
    do j=1,nPhotoBins
       kavg = 0d0
       ksum = 0d0
       do i=1,icount-2
          ! integrate only in the bin range
          if(tmp_energy(i)>=photoBinEleft(j) &
               .and. tmp_energy(i+1)<=photoBinEright(j)) then
             ! numerator integral Jdraine(E)kabs(E)/E
             f1 = tmp_data(i)*Jdraine(i)/tmp_energy(i)
             f2 = tmp_data(i+1)*Jdraine(i+1)/tmp_energy(i+1)
             kavg = kavg + (f1+f2) / 2d0 &
                  * (tmp_energy(i+1)-tmp_energy(i))

             ! denominator integral Jdraine(E)/E
             f1 = Jdraine(i)/tmp_energy(i)
             f2 = Jdraine(i+1)/tmp_energy(i+1)
             ksum = ksum + (f1+f2) / 2d0 &
                  * (tmp_energy(i+1)-tmp_energy(i))
          end if
!!$           if(tmp_energy(i)<photoBinEmid(j) &
!!$                .and. tmp_energy(i+1)>photoBinEmid(j)) then
!!$              kavg = (photoBinEmid(j) - tmp_energy(i)) &
!!$                   / (tmp_energy(i+1) - tmp_energy(i)) &
!!$                   * (tmp_data(i+1) - tmp_data(i)) + tmp_data(i)
!!$              print *,photoBinEmid(j), kavg
!!$           end if
       end do
       ! ratio of the integral is average absorption in the bin
       find_Av_draine_kabs(j) = kavg / (ksum+1d-40)
    end do

  end subroutine find_Av_load_kabs2

  ! *****************************
  ! load kabs from file
  subroutine find_Av_load_kabs(file_name)
    use krome_commons
    use krome_constants
    implicit none
    character(len=*),intent(in),optional::file_name
    character(len=200)::fname
    real*8::opacityDust_org(nPhotoBins)

    ! check if energy bins are set
    if(maxval(photoBinEleft)==0d0) then
      print *, "ERROR: to load kabs for Av G0 finder you"
      print *, " have to initialize some energy bins!"
      stop
    end if

    ! check if optional argument is present
    fname = "kabs_draine_Rv31.dat"
    if(present(file_name)) then
      fname = trim(file_name)
    end if

    opacityDust_org = opacityDust

    call load_opacity_table(fname, "micron")

    ! ratio of the integral is average absorption in the bin
    find_Av_draine_kabs = opacityDust

    opacityDust = opacityDust_org

  end subroutine find_Av_load_kabs

  ! *********************************
  ! given the current photo bin intensity distribution
  ! estimates G0 and Av using the bins in the Draine range
  subroutine estimate_G0_Av(G0, Av, n, d2g)
    use krome_constants
    use krome_commons
    use krome_getphys
    use krome_subs
    implicit none
    real*8,intent(out)::G0, Av
    real*8,intent(in)::d2g, n(nspec)
    real*8::lnG0, mu, ntot
    real*8::ydata(nPhotoBins)
    integer::i
    logical, save ::first_call=.true.
    real*8,  save ::XH,Jdraine(nPhotoBins),xdata(nPhotoBins)
    integer, save ::lb,ub,ndraine
    !$omp threadprivate(first_call,XH,Jdraine,xdata,lb,ub,ndraine)

    if (first_call) then
       ! get non-attenuated draine flux
       Jdraine(:) = get_draine(photoBinEmid(:))

       ! only consider bins that have non-attenuated draine radiation
       lb=0
       do i=1,nPhotoBins
          if (Jdraine(i)>0 .and. lb==0) lb=i
          if (Jdraine(i)>0) ub=i
       end do
       ndraine = ub - lb + 1

       ! mean molecular weight and gas density
       mu = get_mu(n(:))
       ntot = sum(n(1:nspec))

       ! find mass fraction of H-nuclei as xH = mp * n_H / rho
       xH = (p_mass * get_Hnuclei(n(:))) / (p_mass * mu * ntot)

       ! now we can calculate the xdata, which are constant:

       ! loop on photo bins
       do i=lb,ub
          ! compute x in y = Av*x + ln(G0)
          xdata(i-lb+1) = -find_Av_draine_kabs(i) * 1.8d21 * p_mass * d2g / xH
       end do

       ! make sure we only do this once
       ! NOTICE: we assume hydrogen mass fraction is constant
       ! NOTICE: but this is implicitly the case anyway
       ! NOTICE: because below we translate between
       ! NOTICE: column density and Av
       first_call = .false.
    end if

    ! loop on photo bins
    do i=lb,ub
       ! compute y in y = Av*x + ln(G0)
       ydata(i - lb + 1) = log(photoBinJ(i) + 1d-200) - log(Jdraine(i))
    end do

    ! needs at least one bin
    if(ndraine<=1) then
       print *,"ERROR: you want to estimate G0 and Av with less than 2 bins in the"
       print *," Draine range, 5-13.6 eV! Nbins(Draine)=",ndraine
       stop
    end if

    ! call least squares to compute Av and ln(G0)
    call llsq(ndraine, xdata(1:ndraine), ydata(1:ndraine), &
         Av, lnG0)

    ! Apply prior
    if(lnG0 < -7d0 .or. lnG0 > 7d0) then
      if(lnG0 < -7d0) lnG0 = -7d0
      if(lnG0 > 7d0) lnG0 = 7d0
      Av = sum(( ydata(1:ndraine)-lnG0)*xdata(1:ndraine))/sum(xdata(1:ndraine)**2)
    end if

    if(Av < 0d0) then
      Av = 0d0
      lnG0 = sum( ydata(1:ndraine))/ndraine
    endif

    ! return G0
    G0 = exp(lnG0)

  end subroutine estimate_G0_Av

  ! ************************
  function get_draine(energy_list) result(Jdraine)
    use krome_commons
    use krome_constants
    implicit none
    integer::i
    real*8,intent(in)::energy_list(:)
    real*8::x, Jdraine(size(energy_list))

    do i=1,size(energy_list)
       x = energy_list(i) !eV
       !eV/cm2/sr
       if(x<13.6d0.and.x>5d0) then
          Jdraine(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
               * x *planck_eV
       else
          Jdraine(i) = 0d0
       end if
    end do
  end function get_draine

#ENDIFKROME

end module krome_photo
