!################################################################
!PDR test, similar to Gong, Ostriker & Wolfire 2017
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
  use krome_dust, ONLY : compute_Semenov_Tdust
  implicit none
  integer,parameter::nz=1
  integer,parameter::rstep = 500000
  integer::i,j,ii,ios,jscale,jz,jz2, column_bins, zint, NoColumnBins
  real*8::rhogas,m(krome_nspec)
  real*8::tff,ertol,eatol,max_time,t_tot
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea),kkk(krome_nspec),ColumnTot,ColumnTotMax,ColumnTotMin,ColumnLast,dColumn,ColumnFactor
  real*8::Av,heats(krome_nheats),crate,crate_0,NH_cum,NH2_cum,NC_cum, NCO_cum
  real*8::ionH,dissH2,ionC,dissCO,chiFUV,chiLW,chiPE,chi0,dustHeatingRate
  logical::stop_next, converged
  character(len=20) :: filename, zint_str
  real*8, parameter :: Lshield_0 = 1.5428402399039558e+19, a = 0.7, n_0 = 100.0, sigmaD_LW = 1.5e-21, sigmaD_PE = 0.86e-21, bfive=1d0
  real*8 :: Lshield, Nshield, t_cool, ntot_val, ntotchange_cum
  real*8, parameter :: J_FUV_ISRF = 2.1e-4, dustUV_crossSection = 1.e-21

  !zs = (/1d-6, 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 1d0/) !list of metallicities relative to solar
  zs = (/1d0/)

  !set the scaled FUV intensity
  chi0 = 1d0
  !Set the cosmic ray rate, proportional to the FUV intensity; default for ISRF 2x10^-16 s^-1
  crate_0 = 1.e-16 !Note: this is the primary+secondary CR ionization of H2 ~ 0.5 times the primary H rate


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
    filename = trim('PDR_Z') // trim(zint_str)
    filename = trim(filename)
    !Open file
    open(unit=22,file=filename,status='replace',action='write')
    write(22, '(A)', ADVANCE='NO') "#ntot rho Tgas Tdust ColumnTot"
    write(22, '(A)') trim(krome_get_names_header())
    
    print *, 'Metallicity: ', zs(jz2), ' of Solar'

    !INITIAL CONDITIONS
    krome_redshift = 0d0    !redshift
    Tgas = 3d2             !temperature, K
    ntot = 1d2             ! Fixed density of 100cm^-3
    ColumnTotMin = 1d17   ! Minimum column density    
    ColumnTotMax = 1d22   ! Maximum column density
    NoColumnBins = 500   ! Number of column bins in log space
    ColumnFactor = 10**((log10(ColumnTotMax) - log10(ColumnTotMin))/NoColumnBins) ! Column factor in log space
    NH_cum = 0d0           ! Cumulative H column (for shielding)
    NH2_cum = 0d0          ! Cumulative H2 column (for shielding)
    NC_cum = 0d0           ! Cumulative C column (for shielding)
    NCO_cum = 0d0          ! Cumulative CO column (for shielding)

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
    ertol = 1d-4  ! relative min change in a species

    ColumnTot = ColumnTotMin

    do column_bins = 1, 10000

      !species default, cm-3
      x(:) = 1d-40

      !set individual species
      x(KROME_idx_H)         = ntot - 2*1d-6*ntot - 1d-4*ntot
      x(KROME_idx_H2)        = 1d-6*ntot
      x(KROME_idx_E)         = 1d-4*ntot
      x(KROME_idx_Hj)        = 1d-4*ntot
      x(KROME_idx_HE)        = 0.0775*ntot
      x(KROME_idx_Cj)        = 1.6d-4*zs(jz2)*ntot !C is fully ionized
      x(KROME_idx_O)         = 3.2d-4*zs(jz2)*ntot !O is fully neutral
      x(KROME_idx_D)         = 3d-5*ntot

      call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

      !set H2 dissociation reaction rate coeff
      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      ni(:) = n(:)

      !Set shielded quantities and rates
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Nshield = NH_cum + 2*NH2_cum
      Av = Nshield * zs(jz2) / 1.87d21
      call krome_set_user_Av(Av)

      !set H ionization reaction rate coeff
      ionH = 0.0 !No H ionization
      call krome_set_user_ionH(ionH)

      !LW and PE rates
      chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
      chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
      !Dissociation rates
      dissH2 = 5.60d-11*chiLW*get_fshield_H2(NH2_cum,bfive)  !H2 dissociation rate accounting for self-shielding
      call krome_set_user_dissH2(dissH2)
      ionC = 3.1d-10*krome_get_user_is_metal()*chiLW*get_fshield_C(NH2_cum,NC_cum,bfive)
      dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW*get_fshield_CO(NH2_cum,NCO_cum,bfive)
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

      ntotchange_cum = 0d0

      !loop on density steps
      do i=1, rstep

        if (i .eq. 1) then
          !free-fall time, s
          tff = krome_get_free_fall_time(x(:))
          user_tff = tff         !store user tff. NEED THIS LINE FOR SOME WEIRD REASON
        endif

        !print *, "Total number density in new loop: ", sum(x(:))

        !if you do not conserve electrons, the electron abundance will soon go to 0.00
        ntot_val = sum(x(:))
        x(krome_idx_e) = krome_get_electrons(x(:))
        !Now do a renormalization of the number densities since the above changes number densities
        do j=1,krome_nspec
          x(j) = x(j) * ntot_val / sum(x(:))
        end do

        !Set shielded quantities and rates
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Nshield = NH_cum + 2*NH2_cum
        Av = Nshield * zs(jz2) / 1.87d21
        call krome_set_user_Av(Av)

        !set H ionization reaction rate coeff
        ionH = 0.0 !No H ionization
        call krome_set_user_ionH(ionH)

        !LW and PE rates
        chiLW = chi0 * exp(-sigmaD_LW * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
        chiPE = chi0 * exp(-sigmaD_PE * zs(jz2) * Nshield) !Dust extinction, where D linearly scales with Z
        !Dissociation rates
        dissH2 = 5.60d-11*chiLW*get_fshield_H2(NH2_cum,bfive)  !H2 dissociation rate accounting for self-shielding
        call krome_set_user_dissH2(dissH2)
        ionC = 3.1d-10*krome_get_user_is_metal()*chiLW*get_fshield_C(NH2_cum,NC_cum,bfive)
        dissCO = 2.592d-10*krome_get_user_is_metal()*chiLW*get_fshield_CO(NH2_cum,NCO_cum,bfive)
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
        !Absorption rate of UV photons by dust (erg s^-1)
        dustHeatingRate = chiFUV*J_FUV_ISRF*4*pi*dustUV_crossSection*zs(jz2)
        call krome_set_dustheatRad(dustHeatingRate)
        call compute_Semenov_Tdust(x(:), Tgas)
        Tdust = krome_get_Semenov_Tdust()

        ni(krome_idx_Tgas) = Tgas

        !call equilibrium solver
        call krome_equilibrium_xT(x(:),Tgas,dt)

        !avoid negative species
        do ii=1,krome_nmols
          n(ii) = max(x(ii),0d0)
        end do
        n(krome_idx_Tgas) = Tgas

        ! check if we have converged by comparing the error in any species with an relative abundance above eatol
        ! print *, 'haha: ', n(krome_idx_Tgas), ni(krome_idx_Tgas), abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas)
        ! converged = abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .le. ertol &
        !              .or. t_tot .gt. max_time

        ! check if we have converged by comparing the error in any species with an relative abundance above eatol
        converged = maxval(abs(n(1:krome_nmols) - ni(1:krome_nmols)) / max(n(1:krome_nmols),eatol*sum(n(1:krome_nmols)))) .lt. ertol &
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
                    "CONVERGED; Column = ", ColumnTot, " Tgas = ", Tgas, " t_tot/Myr = ", &
                    t_tot/(seconds_per_year*1.e6), " dt = ", dt/(seconds_per_year*1.e6), &
                    " t_cool = ", t_cool/(seconds_per_year*1.e6)
          exit
        endif
      end do

      if(t_tot > max_time .or. abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .gt. ertol) then
        print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years. Reldiff: ', abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas)
        print *, 'Tgas :', n(krome_idx_Tgas)
      end if

      m = get_mass()
      rhogas = sum(n(1:krome_nmols)*m(1:krome_nmols))
      write(22,'(99E17.8e3)') sum(n(1:krome_nmols)),rhogas,Tgas,Tdust,ColumnTot,n(1:krome_nmols)/sum(n(1:krome_nmols))
      flush(22)

      if (stop_next) exit

      !Store last column
      ColumnLast = ColumnTot
      !increase column by the appropriate factor for the next bin
      ColumnTot = ColumnTot * ColumnFactor
      !Change in column
      dColumn = ColumnTot - ColumnLast
      !Add to cumulative columns of H, H2, C, CO
      NH_cum = NH_cum + (n(KROME_idx_H)/sum(n(1:krome_nmols))) * dColumn
      NH2_cum = NH2_cum + (n(KROME_idx_H2)/sum(n(1:krome_nmols))) * dColumn
      NC_cum = NC_cum + (n(KROME_idx_Cj)/sum(n(1:krome_nmols))) * dColumn
      NCO_cum = NCO_cum + (n(KROME_idx_CO)/sum(n(1:krome_nmols))) * dColumn
      !break when max density reached
      if (ColumnTot .gt. ColumnTotMax) then
        ColumnTot = ColumnTotMax
        stop_next = .true.
      endif
    end do

    !write(22, *)

    !Close files
    close(22)
  end do

  !call cool_DustGRREC(5d0,1.e3)

  !say goodbye
  print *,"To plot in python:"
  print *,"ipython> run plot.py"
  print *,"That's all! have a nice day!"

contains

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
    z(:, 1) = (/ 1d0, 8.080d-1, 5.250d-1, 2.434d-1, 5.467d-2, 1.362d-2, 3.378d-3, 5.240d-4 /)
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

end program test_krome_eqbm

