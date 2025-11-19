!################################################################
!Same as the popsicle_semenov_photo_cr test but with 
!GOW network at metallicities >= 0.01*Solar.
!For additional details, see Gong et al. 2017, Hunter et al. 2023
!Kim et al. 2023 ApJS, and the KROME paper (Grassi et al. 2014).
!Author: Piyush Sharda (Leiden, 2025)
!Email: sharda@strw.leidenuniv.nl
!################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_cooling
  use krome_heating
  use krome_getphys
  use krome_phfuncs
  use krome_dust, ONLY : compute_Semenov_Tdust
  use krome_constants
  implicit none
  integer,parameter::nz=1
  integer,parameter::rstep = 500000
  integer::i,unit,ios,jscale,jz,jz2,tendbytff,ii
  real*8::dtH,deldd,rhogas,m(krome_nspec)
  real*8::tff,dd,dd1,ertol,eatol
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools),time, t_cool
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea)
  real*8::Av,heats(krome_nheats),crate,NH,NHj,NH2,chi0, sum_xi, sum_x, Hnuclei_i, Hnuclei
  real*8::ionH,dissH2,ionC,chiCO,chiFUV,chiLW,chiPE,dustHeatingRate,dissCO,Nshield,Lshield,crate_0
  logical::uv_attenuation, cr_attenuation, converged
  real*8, parameter :: J_FUV_ISRF = 2.1e-4, dustUV_crossSection = 1.e-21
  real, parameter :: Lshield_0 = 1.5428402399039558e+19, a = 0.7, n_0 = 100.0, sigmaD_LW = 1.5e-21, sigmaD_PE = 0.86e-21

  zs = (/1d0/) !list of metallicities relative to solar

  !INITIAL CONDITIONS
  krome_redshift = 0d0    !redshift
  ntot = 1d2            !total number density, cm-3
  Tgas = 1d4              !temperature, K
  tendbytff = 100             !end time, s
  !set to True to switch on cosmic ray/UV attenuation
  uv_attenuation = .true.
  cr_attenuation = .true.

  !set chiFUV for photoreactions
  chi0 = 1d0
  crate_0 = 2d-16

  !output header
  write(22, '(A)', ADVANCE='NO') "#time ntot rhotot Tgas"
  write(22, '(A)') trim(krome_get_names_header())

  write(31, '(A)', ADVANCE='NO') "#ntot Tgas sum(cools)"
  write(31, '(A)') trim(krome_get_cooling_names_header())

  write(911, '(A)', ADVANCE='NO') "#ntot Tgas sum(heats)"
  write(911, '(A)') trim(krome_get_heating_names_header())

  !loop over size(zs)*2 so that every second loop is skipped, so that an empty line is created in the output fort.22 file
  !this line break in the output file can then be used to read in output for each zs separately
  do jz = 1,size(zs)*2

     jz2 = (jz+1)/2
     jscale = mod(jz,2)
     if(jscale==0) cycle

    print *, 'Metallicity: ', zs(jz2), ' of Solar'

    call krome_set_zredshift(krome_redshift)
    call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
    call krome_set_metallicity(zs(jz2))
    call krome_set_dust_to_gas(zs(jz2))

    if (zs(jz2) > 0d0) then
      !turn on photo/cr reactions that include metals
      call krome_set_user_is_metal(1d0)
    else
      !turn off photo/cr reactions that include metals
      call krome_set_user_is_metal(0d0)
    endif

    !initialize KROME (mandatory)
    call krome_init()

    !species default, cm-3
    x(:) = 1d-40

    !set individual species
    !Fully atomic initial condition
    x(KROME_idx_H)         = ntot* (1d0 - (2*1d-3 + 3*2.681411d-07 + 1d-4))
    x(KROME_idx_H2)        = 2*1d-3*ntot
    !Fully molecular initial condition
    !x(KROME_idx_H)         = 1d-4*ntot
    !x(KROME_idx_H2)        = 0.5*ntot - ntot*(3*2.681411d-07 + 1d-4)
    x(KROME_idx_E)         = 1.6d-4*zs(jz2)*ntot + 1d-4*ntot + 2.681411e-07*ntot
    x(KROME_idx_Hj)        = 1d-4*ntot
    x(KROME_idx_HE)        = 0.1*ntot
    x(KROME_idx_Cj)        = 1.6d-4*zs(jz2)*ntot !C is fully ionized
    x(KROME_idx_O)         = 3.2d-4*zs(jz2)*ntot !O is fully neutral
    x(KROME_idx_D)         = 3d-5*ntot
    x(KROME_idx_H3j)       = 3*2.681411e-07*ntot

    call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

    !initial Hnuclei
    Hnuclei_i = get_Hnuclei(x(:))

    !set initial density
    dd = ntot

    !open file to write explore data
    open(newunit=unit,file="explore.dat",status="replace")

    print *,"solving..."
    print '(a8,8a11)',"step","time (s)","n(cm-3)","Tgas(K)", "crate", "Av","chiLW", "chiFUV"

    time = 0d0

    ertol = 1d-8  ! relative min change in a species
    eatol = 1d-10 ! absolute min change in a species
    t_cool = 1d10

    !loop on density steps
    do i = 1,rstep

       ni(:) = n(:)

       !store old density
       dd1 = dd

       !free-fall time, s
       tff = krome_get_free_fall_time(x(:))
       user_tff = tff         !store user tff
       dtH = min(0.01d0 * tff,t_cool)     !define time-step
       time = time + dtH

       !if you do not conserve electrons, the electron abundance will soon go to 0.00
       sum_xi = sum(x(1:krome_nmols))
       x(krome_idx_e) = krome_get_electrons(x(:))
       sum_x = sum(x(1:krome_nmols))
       x(1:krome_nmols) = x(1:krome_nmols) * sum_xi / sum_x

       !set time-step
       dt = dtH
       !initial Hnuclei
       Hnuclei = get_Hnuclei(x(:))
       !Shielding
       if(uv_attenuation .or. cr_attenuation) then
        Lshield = Lshield_0 * (Hnuclei/n_0)**(-a)
        Nshield = Lshield * Hnuclei
       else
        Nshield = 0.0
       endif
       Av = Nshield * zs(jz2) / 1.87d21
       call krome_set_user_Av(Av)

       if(cr_attenuation) then
        if(Nshield < 9.35e20) then
          crate = crate_0
        else
          crate = crate_0 * (Nshield/9.35e20)**(-1)
        endif
       else
        crate = crate_0
       endif

       call krome_set_user_crate(crate)

       ionH = 0.0 !No H ionization
       call krome_set_user_ionH(ionH)
       !LW and PE rates
       chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
       chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
       !Dissociation rates
       dissH2 = 5.60d-11*chiLW
       ionC = 3.1d-10*krome_get_user_is_metal()*chiLW
       dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW

       if(uv_attenuation) then
        dissH2 = dissH2 * get_fshield_H2(Nshield * x(KROME_idx_H2)/Hnuclei, 1d0)
        ionC = ionC * get_fshield_C(Nshield * x(KROME_idx_H2),Nshield * x(KROME_idx_Cj)/Hnuclei, 1d0)
        dissCO = dissCO * get_fshield_CO(Nshield * x(KROME_idx_H2),Nshield * x(KROME_idx_CO)/Hnuclei, 1d0)
       end if

       call krome_set_user_dissH2(dissH2)
       call krome_set_user_ionC(ionC)
       call krome_set_user_dissCO(dissCO)

       ! Case of TIGRESS NCR (comment above two lines and uncomment below)
       !chiCO = krome_get_user_is_metal()*chiLW
       !call krome_set_user_chiCO(chiCO)


       !FUV rate for photoelectric heating (FUV = LW + PE; both of these are attenuated separately as above)
       chiFUV = (chiPE * 1.8e-4 + chiLW * 3.e-5)/(2.1e-4) !Scale and sum attenuated ISRF LW/PE intensities to the mean FUV intensity
       call krome_set_chiFUV(chiFUV)

       !Absorption rate of UV photons by dust (erg s^-1)
       dustHeatingRate = chiFUV*J_FUV_ISRF*4*pi*dustUV_crossSection*zs(jz2)
       call krome_set_dustheatRad(dustHeatingRate)
       call compute_Semenov_Tdust(x(:), Tgas)
       Tdust = krome_get_Semenov_Tdust()

       !dump cooling rates for Tgas going into the calculation
       n(1:krome_nmols) = x(:)
       n(KROME_idx_Tgas) = Tgas
       cools(:) = get_cooling_array(n(:),Tgas)
       write(31,'(99E14.5e3)') time, dd,Tgas, sum(cools), cools(:)
       kk(:) = krome_get_coef(Tgas,x(:))
       heats(:) = get_heating_array(n(:),Tgas,kk(:),0d0) !TODO: pass nH2dust instead of 0d0 as the third argument
       write(911,'(99E14.5e3)') time, dd,Tgas, sum(heats), heats(:)

       !returns to user array
       x(:) = n(1:krome_nmols)

       ni(krome_idx_Tgas) = Tgas

       !solve the chemistry
       call krome(x(:),Tgas,dt)

       !avoid negative species
        do ii=1,krome_nmols
          n(ii) = max(x(ii),0d0)
        end do

       !Rescale number densities to conserve H nuclei
       Hnuclei = get_Hnuclei(n(:))
       n(1:krome_nmols) = n(1:krome_nmols) * Hnuclei_i/Hnuclei
       Hnuclei = get_Hnuclei(n(:))
       n(krome_idx_Tgas) = Tgas

       converged = maxval(abs(n(1:krome_nmols) - ni(1:krome_nmols)) / max(n(1:krome_nmols),eatol*sum(n(1:krome_nmols)))) .lt. ertol &
          .or. time .gt. tendbytff * user_tff

       t_cool = (Hnuclei * boltzmann_erg * Tgas)/(cooling(n(:),Tgas))


       !print some output
       m = get_mass()
       rhogas = sum(x(:)*m(1:krome_nmols))
       write(22,'(99E17.8e3)') time,sum(x(:)),rhogas,Tgas,x(:)
       if(mod(i,100)==0) then
          !totheat = krome_get_heating(x(:), Tgas)
          !totcool = krome_get_heating(x(:), Tgas)
          print '(I6,30E11.3)',i,time,dd,Tgas,crate,krome_get_user_Av()
          call krome_print_best_flux(x(:),Tgas,5)
          call krome_explore_flux(x(:),Tgas,unit,dd)
       end if

       if(converged) exit
    end do
    write(22,*)
    write(31,*)
    write(911,*)
  end do

  !close explore data file
  close(unit)

  !say goodbye
  print *,"To plot in python:"
  print *,"ipython> run plot.py"
  print *,"That's all! have a nice day!"

contains

  !CR attenutation following Table F1 of Padovani et al. 2018, A&A 614, A111
  real(8) function calculate_F(b)
    implicit none
    real(8), intent(in) :: b
    real(8) :: log_b
    integer :: k

    ! Coefficients c_k
    real(8), dimension(0:9) :: coefficients
    data coefficients / &
         1.001098610761d7, -4.231294690194d6, 7.921914432011d5, -8.623677095423d4, &
         6.015889127529d3, -2.789238383353d2, 8.595814402406d0, -1.698029737474d-1, &
         1.951179287567d-3, -9.937499546711d-6 /

    ! Calculate log10(b)
    log_b = log10(b)

    ! Initialize F to zero
    calculate_F = 0.0d0

    ! Calculate F using the polynomial equation
    do k = 0, 9
      calculate_F = calculate_F + coefficients(k) * (log_b ** k)
    end do

  end function calculate_F

  !
  !===============================================================================
  !
  function get_fshield_H2(NH2,bfive)

  !
  ! Returns the H2 self-shielding function
  ! Eq 12 of Wolcott-Green, Haiman and Bryan 2011: note this is slightly different from DB function
    implicit none
    real*8, intent(in) :: NH2, bfive
    real*8 :: get_fshield_H2

    get_fshield_H2 = 0.965/(1+(NH2/(5.e14*bfive)))**1.1 + &
                    0.035/((1. + (NH2/5.e14))**0.5) * exp(-8.5 * 1.e-4 * (1. + (NH2/5.e14))**0.5)
    return
  end function get_fshield_H2

  function get_fshield_H(NH,bfive)
  !
  ! Returns the self-shielding due to Lyman-alpha lines on the LW band
  ! Eq 15 of Wolcott-Green, Haiman and Bryan 2011
    implicit none
    real*8, intent(in) :: NH, bfive
    real*8 :: get_fshield_H

    get_fshield_H = max( 1/(1+(NH/(2.85e23)))**1.6 * exp(-0.15 * (NH/(2.85e23))), 1e-15 )
    return
  end function get_fshield_H

  function get_fshield_C(NH2,NC,bfive)
  !
  ! Returns the shielding factor for C using the treatment of Tielens & Hollenbach (1985)
  ! Eq 9 in Gong, Ostriker & Wolfire 2017
    implicit none
    real*8, intent(in) :: NH2, NC, bfive
    real*8 :: get_fshield_C

    get_fshield_C = exp(-NC * 1.6e-17) * exp(-NH2 * 2.8e-22)/(1 + (2.8e-22 * NH2))
    return
  end function get_fshield_C

  !*****************************
    !2D interpolation at (x0,y0) for (x(:), y(:)) in z(:,:)
    !Added by Piyush Sharda in 2024 for CO shielding
  function interpolate2DCO(x, y, z, x0, y0)
    implicit none
    real*8 :: interpolate2DCO
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
    interpolate2DCO = (1 - t) * (1 - u) * z(i, j) + t * (1 - u) * z(i + 1, j) + &
        (1 - t) * u * z(i, j + 1) + t * u * z(i + 1, j + 1)

  end function interpolate2DCO

  function get_fshield_CO(NH2,NCO,bfive)
  !
  ! Returns the shielding factor for CO using tabulated data from Visser et al. 2009, compiled by Gong et al. 2017
  ! Tabulated data procured from: https://github.com/munan/pdr/blob/master/shielding.cpp
    implicit none
    real*8, intent(in) :: NH2, NCO, bfive
    real*8 :: get_fshield_CO
    real*8 :: x(8), y(6), z(8, 6), clipped_x,clipped_y

    x = (/ 0d0, 13d0, 14d0, 15d0, 16d0, 17d0, 18d0, 19d0 /) !N_CO
    y = (/ 0d0, 19d0, 20d0, 21d0, 22d0, 23d0 /) !N_H2
    z(:, 1) = (/ 1d0, 8.080d-1, 5.250d-1, 2.434d-1, 5.467d-2, 1.362d-2, 3.378d-3, 5.240d-5 /)
    z(:, 2) = (/ 8.176d-1, 6.347d-1, 3.891d-1, 1.787d-1, 4.297d-2, 1.152d-2, 2.922d-3, 4.662d-4 /)
    z(:, 3) = (/ 7.223d-1, 5.624d-1, 3.434d-1, 1.540d-1, 3.515d-2, 9.231d-3, 2.388d-3, 3.899d-4 /)
    z(:, 4) = (/ 3.260d-1, 2.810d-1, 1.953d-1, 8.726d-2, 1.907d-2, 4.768d-3, 1.150d-3, 1.941d-4 /)
    z(:, 5) = (/ 1.108d-2, 1.081d-2, 9.033d-3, 4.441d-3, 1.102d-3, 2.644d-4, 7.329d-5, 1.437d-5 /)
    z(:, 6) = (/ 3.938d-7, 3.938d-7, 3.936d-7, 3.923d-7, 3.901d-7, 3.893d-7, 3.890d-7, 3.875d-7 /)

    !Clip logNCO and logNH2 to the ranges in the data
    clipped_x = max(x(1), min(log10(NCO), x(8)))
    clipped_y = max(y(1), min(log10(NH2), y(6)))
    get_fshield_CO = 1d1**interpolate2DCO(x, y, log10(z), clipped_x, clipped_y)
    return
  end function get_fshield_CO


end program test_krome
