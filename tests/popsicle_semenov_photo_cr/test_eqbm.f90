!################################################################
!Same as the popsicle_semenov_photo_cr test but keeps density constant
!and evolves the temperature to equilibrium
!For additional details, see Omukai 2000, +2005,
!Glover et al. 2010, Kim et al. 2023 ApJS, and
!the KROME paper (Grassi et al. 2014).
!Author: Piyush Sharda (Leiden, 2024)
!Email: sharda@strw.leidenuniv.nl
!################################################################
program test_krome_eqbm

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_cooling
  use krome_heating
  use krome_getphys
  use krome_phfuncs
  use krome_constants
  implicit none
  integer,parameter::nz=1
  integer,parameter::rstep = 500000
  integer::i,ii,ios,jscale,jz,jz2, dens_bins
  real*8::rhogas,m(krome_nspec)
  real*8::tff,ertol,eatol,max_time,t_tot
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea),kkk(krome_nspec)
  real*8::Av,heats(krome_nheats),crate,NH,NHj,NH2
  real*8::ionH,dissH2,ionC,dissCO,chiFUV
  logical::crate_attenuation, stop_next, converged

  !zs = (/1d-6, 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 1d0/) !list of metallicities relative to solar
  zs = (/1d0/)

  !set to True to switch on cosmic ray attenuation
  crate_attenuation = .False.

  !set chiFUV for photoreactions
  chiFUV = 1d0


  ! Switches to decide when equilibrium has been reached
  ertol = 1d-8  ! relative min change in a species
  max_time=seconds_per_year*5d8 ! max time we will be integrating for

  !output header
  write(22, '(A)', ADVANCE='NO') "#ntot rhotot Tgas Tdust"
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
    Tgas = 3d2              !temperature, K
    ntot = 1d-2

    call krome_set_zredshift(krome_redshift)
    call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
    call krome_set_metallicity(zs(jz2))
    call krome_set_chiFUV(chiFUV)

    if (zs(jz2) > 0d0) then
      !turn on photo/cr reactions that include metals
      call krome_set_user_is_metal(1d0)
    else
      !turn off photo/cr reactions that include metals
      call krome_set_user_is_metal(0d0)
    endif

    !initialize KROME (mandatory)
    call krome_init()

    !Cosmic ray ionization rate
    if (crate_attenuation) then
      !Attenuate following Appendix F of Padovani et al. 2018
      crate = 10**calculate_F(NH + NHj + 2d0*NH2)
      print *, 'Using cosmic ray attenutation'
    else
      crate = 2d-16
      print *, 'Using constant cosmic ray ionization rate'
    endif
    print *, 'Initial crate: ', crate
    call krome_set_user_crate(crate)

    !switch to tell when to stop the calculation
    stop_next = .false.

    do dens_bins = 1, 10000

      !species default, cm-3
      x(:) = 1d-40

      !set individual species
      x(KROME_idx_H)         = ntot
      x(KROME_idx_H2)        = 1d-6*ntot
      x(KROME_idx_E)         = 1d-4*ntot
      x(KROME_idx_Hj)        = 1d-4*ntot
      x(KROME_idx_HE)        = 0.0775*ntot
      x(KROME_idx_Cj)        = 0.927d-4*zs(jz2)*ntot !C is fully ionized
      x(KROME_idx_O)         = 3.568d-4*zs(jz2)*ntot !O is fully neutral

      call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

      !set initial Av, following equation 3 of Gong, Ostriker and Wolfire 2017
      NH  = krome_num2col(x(KROME_idx_H), x(:), Tgas)
      NHj = krome_num2col(x(KROME_idx_Hj), x(:), Tgas)
      NH2 = krome_num2col(x(KROME_idx_H2), x(:), Tgas)
      Av = (NH + NHj + 2d0*NH2) *zs(jz2) / 1.87d21
      call krome_set_user_Av(Av)
      print *, 'numdens: ', sum(x(:))

      !set H ionization reaction rate coeff
      ionH = 2.19d-12*exp(-1.14e4*Av)*chiFUV
      call krome_set_user_ionH(ionH)
      !set H2 dissociation reaction rate coeff
      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      ni(:) = n(:)

      dissH2 = 5.60d-11*exp(-3.74*Av)*krome_fshield(n,Tgas)*chiFUV
      call krome_set_user_dissH2(dissH2)
      ionC = 3.1d-10*exp(-3.*Av)*krome_fshield_C(n,Tgas)*krome_get_user_is_metal()*chiFUV
      dissCO = 2.d-10*exp(-3.53*Av)*krome_fshield_CO(n,Tgas)*krome_get_user_is_metal()*chiFUV
      call krome_set_user_ionC(ionC)
      call krome_set_user_dissCO(dissCO)

      dt = seconds_per_year * 1d2
      t_tot = dt
      converged = .false.

      !loop on density steps
      do i=1, rstep

        if (i .eq. 1) then
          !free-fall time, s
          tff = krome_get_free_fall_time(x(:))
          user_tff = tff         !store user tff. NEED THIS LINE FOR SOME WEIRD REASON
        endif

        !if you do not conserve electrons, the electron abundance will soon go to 0.00
        x(krome_idx_e) = krome_get_electrons(x(:))

        NH  = krome_num2col(x(KROME_idx_H), x(:), Tgas)
        NHj = krome_num2col(x(KROME_idx_Hj), x(:), Tgas)
        NH2 = krome_num2col(x(KROME_idx_H2), x(:), Tgas)
        Av = (NH + NHj + 2d0*NH2) *zs(jz2)/ 1.87d21
        call krome_set_user_Av(Av)

        if (crate_attenuation) then
          crate = 10**calculate_F(NH + NHj + 2d0*NH2)
          call krome_set_user_crate(crate)
        endif

        !set H ionization reaction rate coeff
        ionH = 2.19d-12*exp(-1.14e4*Av)*chiFUV
        call krome_set_user_ionH(ionH)
        !set H2 dissociation reaction rate coeff
        dissH2 = 5.60d-11*exp(-3.74*Av)*krome_fshield(n,Tgas)*chiFUV
        call krome_set_user_dissH2(dissH2)
        ionC = 3.1d-10*exp(-3.*Av)*krome_fshield_C(n,Tgas)*krome_get_user_is_metal()*chiFUV
        dissCO = 2.d-10*exp(-3.53*Av)*krome_fshield_CO(n,Tgas)*krome_get_user_is_metal()*chiFUV
        call krome_set_user_ionC(ionC)
        call krome_set_user_dissCO(dissCO)

        Tdust = krome_get_Semenov_Tdust()

        !solve the chemistry
        print *, 'krome before: ', sum(x(:)), Tgas
        call krome_equilibrium_xT(x(:),Tgas,dt)
        print *, 'krome after: ', sum(x(:)), Tgas

        !avoid negative species
        do ii=1,krome_nmols
          n(ii) = max(x(ii),0d0)
        end do
        n(krome_idx_Tgas) = Tgas

        kkk(1:krome_nmols) = krome_conserve(n(1:krome_nmols),ni(1:krome_nmols))
        n(:) = kkk(:)

        ! check if we have converged by comparing the error in any species with an relative abundance above eatol
        converged = abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .lt. ertol &
                     .or. t_tot .gt. max_time

        ! Increase integration time by a reasonable factor
        if(.not. converged) then
          dt = dt * 3.
          t_tot = t_tot + dt
          ni = n
        else
          exit
        endif
      end do

      !returns to user array
      x(:) = n(1:krome_nmols)

      if(t_tot > max_time .or. abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) > 0.1) then
        print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
        print *, 'Tgas :', n(krome_idx_Tgas)
      end if
      print *, ' '

      rhogas = sum(x(:)*m(1:krome_nmols))
      write(22,'(99E17.8e3)') sum(x(:)),Tgas,x(:)

      if (stop_next) exit

      !increase density by 2x for the next bin
      ntot = ntot * 2
      !break when max density reached
      if (ntot .gt. 1d6) then
        ntot = 1d6
        stop_next = .true.
      endif
    end do

    write(22, *)
  end do

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

end program test_krome_eqbm

