!################################################################
!Same as the one-zone ISM equilibrium test, but with a prescribed shielding column density
!For additional details, see Sec. 7.2.2 of Kim+23,
!Author: Shyam Menon (CCA/Rutgers, 2024)
!Email: smenon@flatironinstitute.org
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
  integer::i,ii,ios,jscale,jz,jz2, dens_bins, zint
  real*8::rhogas,m(krome_nspec)
  real*8::tff,ertol,eatol,max_time,t_tot
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea),kkk(krome_nspec)
  real*8::Av,heats(krome_nheats),crate,crate_0,NH,NHj,NH2
  real*8::ionH,dissH2,ionC,dissCO,chiFUV,chiLW,chiPE,chi0
  logical::stop_next, converged
  character(len=20) :: filename, zint_str
  real, parameter :: Lshield_0 = 1.5428402399039558e+19, a = 0.7, n_0 = 100.0, sigmaD_LW = 1.5e-21, sigmaD_PE = 0.86e-21
  real :: Lshield, Nshield, t_cool

  !zs = (/1d-6, 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 1d0/) !list of metallicities relative to solar
  zs = (/1d0/)

  !set the scaled FUV intensity
  chi0 = 1d0
  !Set the cosmic ray rate, proportional to the FUV intensity; default for ISRF 2x10^-16 s^-1
  crate_0 = 2d-16 * chi0


  max_time=seconds_per_year*1.e9 ! max time we will be integrating for = 1000 Myrs (1Gyr)

  !loop over size(zs)*2 so that every second loop is skipped, so that an empty line is created in the output fort.22 file
  !this line break in the output file can then be used to read in output for each zs separately
  do jz = 1,size(zs)*2

    jz2 = (jz+1)/2
    jscale = mod(jz,2)
    if(jscale==0) cycle

    !Deduce filename from metallicity
    zint = int(log10(zs(jz2)))
    if(zint == 0) then
      write(zint_str, '(I1)') zint
    else
      write(zint_str, '(I2)') zint
    endif
    filename = trim('AB_Z') // trim(zint_str)
    filename = trim(filename)
    !Open file
    open(unit=22,file=filename,status='replace',action='write')
    write(22, '(A)', ADVANCE='NO') "#ntot rho Tgas Tdust"
    write(22, '(A)') trim(krome_get_names_header())

    filename = trim('COOL_Z') // trim(zint_str)
    filename = trim(filename)
    !Open file
    open(unit=31,file=filename,status='replace',action='write')
    write(31, '(A)', ADVANCE='NO') "#ntot Tgas sum(cools)"
    write(31, '(A)') trim(krome_get_cooling_names_header())

    filename = trim('HEAT_Z') // trim(zint_str)
    filename = trim(filename)
    !Open file
    open(unit=911,file=filename,status='replace',action='write')
    write(911, '(A)', ADVANCE='NO') "#ntot Tgas sum(heats)"
    write(911, '(A)') trim(krome_get_heating_names_header())
    
    print *, 'Metallicity: ', zs(jz2), ' of Solar'

    !INITIAL CONDITIONS
    krome_redshift = 0d0    !redshift
    Tgas = 3d2             !temperature, K
    ntot = 1d-2

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

    print *, 'Initial crate: ', crate_0

    !switch to tell when to stop the calculation
    stop_next = .false.

    ! Switches to decide when equilibrium has been reached
    ertol = 1d-8  ! relative min change in a species

    do dens_bins = 1, 10000

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

      call krome_set_Semenov_Tdust(6d0) !Dust at 6K

      !set H2 dissociation reaction rate coeff
      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      ni(:) = n(:)

      !Set shielded quantities and rates
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Lshield = Lshield_0 * (sum(x(:))/n_0)**(-a)
      Nshield = Lshield * sum(x(:))
      Av = Nshield * zs(jz2) / 1.87d21
      call krome_set_user_Av(Av)

      !set H ionization reaction rate coeff
      ionH = 0.0 !No H ionization
      call krome_set_user_ionH(ionH)

      !LW and PE rates
      chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
      chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
      !Dissociation rates
      dissH2 = 5.60d-11*chiLW*krome_fshield(n,Tgas)
      call krome_set_user_dissH2(dissH2)
      ionC = 3.1d-10*krome_get_user_is_metal()*chiLW*krome_fshield_C(n,Tgas)
      dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW*krome_fshield_CO(n,Tgas)
      call krome_set_user_ionC(ionC)
      call krome_set_user_dissCO(dissCO)

      !FUV rate for photoelectric heating (FUV = LW + PE; both of these are attenuated separately as above)
      chiFUV = (chiPE * 1.8e-4 + chiLW * 3.e-5)/(2.1e-4) !Scale and sum attenuated ISRF LW/PE intensities to the mean FUV intensity
      call krome_set_chiFUV(chiFUV)
      !Shield the CR rate by Eq. 55 of Kim+23
      if(Nshield < 9.35e20) then
        crate = crate_0
      else
        crate = crate_0 * (Nshield/9.35e20)**(-1)
      endif
      call krome_set_user_crate(crate)

      !print *, 'numdens, Av, chiLW, chiPE, chiFUV, crate : ', ntot, Av, chiLW, chiPE, chiFUV, Nshield, crate

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Shielding done

      dt = seconds_per_year * 1d4 !0.1 Myr initial time step
      t_tot = dt
      converged = .false.

      !Higher densities, lower tolerance for convergence
      if(ntot .gt. 1.e2) ertol = 1d-8

      !loop on density steps
      do i=1, rstep

        if (i .eq. 1) then
          !free-fall time, s
          tff = krome_get_free_fall_time(x(:))
          user_tff = tff         !store user tff. NEED THIS LINE FOR SOME WEIRD REASON
        endif

        !if you do not conserve electrons, the electron abundance will soon go to 0.00
        x(krome_idx_e) = krome_get_electrons(x(:))

        !Set shielded quantities and rates
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Lshield = Lshield_0 * (sum(x(:))/n_0)**(-a)
        Nshield = Lshield * sum(x(:))
        Av = Nshield * zs(jz2) / 1.87d21
        call krome_set_user_Av(Av)

        !set H ionization reaction rate coeff
        ionH = 0.0 !No H ionization
        call krome_set_user_ionH(ionH)

        !LW and PE rates
        chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
        chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
        !Dissociation rates
        dissH2 = 5.60d-11*chiLW*krome_fshield(n,Tgas)
        call krome_set_user_dissH2(dissH2)
        ionC = 3.1d-10*krome_get_user_is_metal()*chiLW*krome_fshield_C(n,Tgas)
        dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW*krome_fshield_CO(n,Tgas)
        call krome_set_user_ionC(ionC)
        call krome_set_user_dissCO(dissCO)

        !FUV rate for photoelectric heating (FUV = LW + PE; both of these are attenuated separately as above)
        chiFUV = (chiPE * 1.8e-4 + chiLW * 3.e-5)/2.1e-4 !Scale and sum attenuated ISRF LW/PE intensities to the mean FUV intensity
        call krome_set_chiFUV(chiFUV)
        !Shield the CR rate by Eq. 55 of Kim+23
        if(Nshield < 9.35e20) then
          crate = crate_0
        else
          crate = crate_0 * (Nshield/9.35e20)**(-1)
        endif
        call krome_set_user_crate(crate)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Shielding done
        Tdust = krome_get_Semenov_Tdust()

        ni(krome_idx_Tgas) = Tgas

        !solve the chemistry
        call krome_equilibrium_xT(x(:),Tgas,dt)
        !print *, 'krome after: ', sum(x(:)), Tgas, abs(ni(krome_idx_Tgas) - Tgas) / ni(krome_idx_Tgas)

        !avoid negative species
        do ii=1,krome_nmols
          n(ii) = max(x(ii),0d0)
        end do
        n(krome_idx_Tgas) = Tgas

        kkk(1:krome_nmols) = krome_conserve(n(1:krome_nmols),ni(1:krome_nmols))
        n(:) = kkk(:)
        n(krome_idx_Tgas) = Tgas

        ! check if we have converged by comparing the error in any species with an relative abundance above eatol
        ! print *, 'haha: ', n(krome_idx_Tgas), ni(krome_idx_Tgas), abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas)
        converged = abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .le. ertol &
                     .or. t_tot .gt. max_time

        !Compute cooling time; t_cool = nk_BT/Lambda; where Lambda is in erg cm^-3 s^-1
        t_cool = (sum(x(:)) * boltzmann_erg * Tgas)/(cooling(n(:),Tgas))

        ! Increase integration time by a reasonable factor
        if(.not. converged) then
          !dt = dt * 3.
          dt = MIN(t_cool,dt*3.0)
          t_tot = t_tot + dt
          ni = n
        else
          write (*, '(A, E12.4, A, E12.4, A, E12.4, A, E12.4, A, E12.4)') &
                    "CONVERGED; ntot = ", sum(x(:)), " Tgas = ", Tgas, " t_tot/Myr = ", &
                    t_tot/(seconds_per_year*1.e6), " dt = ", dt/(seconds_per_year*1.e6), &
                    " t_cool = ", t_cool/(seconds_per_year*1.e6)
          exit
        endif
      end do

      !dump cooling rates for Tgas going into the calculation
      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      cools(:) = get_cooling_array(n(:),Tgas)
      write(31,'(99E14.5e3)') ntot, Tgas, sum(cools), cools(:)
      kk(:) = krome_get_coef(Tgas,x(:))
      heats(:) = get_heating_array(n(:),Tgas,kk(:),0d0) !TODO: pass nH2dust instead of 0d0 as the third argument
      write(911,'(99E14.5e3)') ntot, Tgas, sum(heats), heats(:)

      !returns to user array
      x(:) = n(1:krome_nmols)

      if(t_tot > max_time .or. abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .gt. ertol) then
        print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years. Reldiff: ', abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas)
        print *, 'Tgas :', n(krome_idx_Tgas)
      end if

      m = get_mass()
      rhogas = sum(x(:)*m(1:krome_nmols))
      write(22,'(99E17.8e3)') sum(x(:)),rhogas,Tgas,Tdust,x(:)/sum(x(:))

      if (stop_next) exit

      !increase density by 2x for the next bin
      ntot = ntot * 1.2
      !break when max density reached
      if (ntot .gt. 1.e6) then
        ntot = 1.e6
        stop_next = .true.
      endif
    end do

    !write(22, *)

    !Close files
    close(22)
    close(31)
    close(911)
  end do

  !call cool_DustGRREC(5d0,1.e3)

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

