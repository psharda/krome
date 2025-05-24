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
  integer,parameter::nz=3
  integer,parameter::rstep = 500000
  integer::i,unit,ios,jscale,jz,jz2,tendbytff
  real*8::dtH,deldd,rhogas,m(krome_nspec)
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),cools(krome_ncools),time
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea)
  real*8::Av,heats(krome_nheats),crate,NH,NHj,NH2,chi0
  real*8::ionH,dissH2,ionC,chiCO,chiFUV,chiLW,chiPE,dustHeatingRate,dissCO,Nshield,Lshield,crate_0
  logical::crate_attenuation
  real*8, parameter :: J_FUV_ISRF = 2.1e-4, dustUV_crossSection = 1.e-21
  real, parameter :: Lshield_0 = 1.5428402399039558e+19, a = 0.7, n_0 = 100.0, sigmaD_LW = 1.5e-21, sigmaD_PE = 0.86e-21

  zs = (/1d-2, 1d-1, 1d0/) !list of metallicities relative to solar

  !set to True to switch on cosmic ray attenuation
  crate_attenuation = .false.

  !set chiFUV for photoreactions
  chi0 = 1d0
  crate_0 = 1d-16 * chi0

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

    !INITIAL CONDITIONS
    krome_redshift = 0d0    !redshift
    ntot = 5d2              !total density, cm-3
    Tgas = 300              !temperature, K
    tendbytff = 100             !end time, s

    call krome_set_zredshift(krome_redshift)
    call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
    call krome_set_metallicity(zs(jz2))

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
    x(KROME_idx_H)         = ntot
    x(KROME_idx_H2)        = 1d-6*ntot
    x(KROME_idx_E)         = 1d-4*ntot
    x(KROME_idx_Hj)        = 1d-4*ntot
    x(KROME_idx_HE)        = 0.0775*ntot
    x(KROME_idx_Cj)        = 1.6d-4*zs(jz2)*ntot !C is fully ionized
    x(KROME_idx_O)         = 3.2d-4*zs(jz2)*ntot !O is fully neutral

    call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

    !set initial density
    dd = ntot

    !open file to write explore data
    open(newunit=unit,file="explore.dat",status="replace")

    print *,"solving..."
    print '(a8,8a11)',"step","time (s)","n(cm-3)","Tgas(K)", "crate", "Av","chiLW", "chiFUV"

    time = 0d0

    !loop on density steps
    do i = 1,rstep

       !store old density
       dd1 = dd

       !free-fall time, s
       tff = krome_get_free_fall_time(x(:))
       user_tff = tff         !store user tff
       dtH = 0.01d0 * tff     !define time-step
       time = time + dtH

       !if you do not conserve electrons, the electron abundance will soon go to 0.00
       x(krome_idx_e) = krome_get_electrons(x(:))

       !set time-step
       dt = dtH
       !Shielding
       Lshield = Lshield_0 * (sum(x(:))/n_0)**(-a)
       Nshield = Lshield * sum(x(:))
       Av = Nshield * zs(jz2) / 1.87d21
       call krome_set_user_Av(Av)

       if(Nshield < 9.35e20) then
        crate = crate_0
       else
        crate = crate_0 * (Nshield/9.35e20)**(-1)
       endif
       call krome_set_user_crate(crate)

       ionH = 0.0 !No H ionization
       call krome_set_user_ionH(ionH)
       !LW and PE rates
       chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
       chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
       !Dissociation rates
       dissH2 = 5.60d-11*chiLW*krome_fshield(n,Tgas)
       call krome_set_user_dissH2(dissH2)
       ionC = 3.1d-10*krome_get_user_is_metal()*chiLW*krome_fshield_C(n,Tgas)
       call krome_set_user_ionC(ionC)
       dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW*krome_fshield_CO(n,Tgas)
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

       !solve the chemistry
       call krome(x(:),Tgas,dt)

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

       if(time > tendbytff * user_tff) exit
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

end program test_krome
