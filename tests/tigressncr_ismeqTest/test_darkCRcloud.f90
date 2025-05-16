!################################################################
!One zone test for dark cloud CR chemistry, to compare with Bialy & Sternberg 2016
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
  integer,parameter::nz=10000
  integer,parameter::rstep = 500000
  integer::i,ii,ios,jscale,jz,jz2, dens_bins, zint, j
  real*8::rhogas,m(krome_nspec)
  real*8::tff,ertol,eatol,max_time,t_tot, ntot_val, Tgasinit
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust,kk(krome_nrea),kkk(krome_nspec), zmetal
  real*8::Av,heats(krome_nheats),crate,crate_0,NH,NHj,NH2
  real*8::ionH,dissH2,ionC,chiCO,chiFUV,chiLW,chiPE,chi0,dustHeatingRate
  logical::stop_next, converged
  character(len=20) :: filename, zint_str
  real, parameter :: Lshield_0 = 1.5428402399039558e+19, a = 0.7, n_0 = 100.0, sigmaD_LW = 1.5e-21, sigmaD_PE = 0.86e-21
  real :: Lshield, Nshield, t_cool
  real*8, parameter :: J_FUV_ISRF = 2.1e-4, dustUV_crossSection = 1.e-21


  !set the scaled FUV intensity
  chi0 = 0d0

  max_time=seconds_per_year*1.e9 ! max time we will be integrating for = 1000 Myrs (1Gyr)

  !INITIAL CONDITIONS
  krome_redshift = 0d0    !redshift
  Tgasinit = 1d2             !temperature, K
  Tgas = Tgasinit
  ntot = 1d3
  !Set the cosmic ray rate, proportional to the FUV intensity; default for ISRF 10^-16 s^-1
  crate_0 = 1d-16 

  zmetal = 1d-3
  !switch to tell when to stop the calculation
  stop_next = .false.
  ! Switches to decide when equilibrium has been reached
  ertol = 1d-4  ! relative min change in a species
  eatol = 1d-12

  filename = trim('DC_Z')
  filename = trim(filename)
  !Open file
  open(unit=22,file=filename,status='replace',action='write')
  write(22, '(A)', ADVANCE='NO') "#ntot rho Tgas Tdust Zmetal"
  write(22, '(A)') trim(krome_get_names_header())

  call krome_set_zredshift(krome_redshift)
  call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
  do jz=1, nz 

    call krome_set_metallicity(zmetal)

    if (zmetal > 0d0) then
      !turn on photo/cr reactions that include metals
      call krome_set_user_is_metal(1d0)
    else
      !turn off photo/cr reactions that include metals
      call krome_set_user_is_metal(0d0)
    endif

    !initialize KROME (mandatory)
    call krome_init()

    
    print *, 'Metallicity: ', zmetal, ' of Solar'

    !species default, cm-3
    x(:) = 1d-40

    !set individual species
    x(KROME_idx_H)         = ntot - 2*1d-6*ntot - 1d-4*ntot
    x(KROME_idx_H2)        = 1d-6*ntot
    x(KROME_idx_E)         = 1d-4*ntot
    x(KROME_idx_Hj)        = 1d-4*ntot
    !x(KROME_idx_HE)        = 0.0775*ntot
    x(KROME_idx_Cj)        = 1.6d-4*zmetal*ntot !C is fully ionized
    x(KROME_idx_O)         = 3.2d-4*zmetal*ntot !O is fully neutral
    !x(KROME_idx_SIj)       = 1.7d-6*zmetal*ntot !Si is fully ionized

    call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

    n(1:krome_nmols) = x(:)
    n(KROME_idx_Tgas) = Tgas
    ni(:) = n(:)

    !set H ionization reaction rate coeff
    ionH = 0.0 !No H ionization
    call krome_set_user_ionH(ionH)

    !Dissociation rates
    dissH2 = 0d0
    ionC = 0d0
    chiCO = 0d0
    chiFUV = 0d0
    Av = 0.0
    call krome_set_user_dissH2(dissH2)
    call krome_set_user_ionC(ionC)
    call krome_set_user_chiCO(chiCO)
    call krome_set_chiFUV(chiFUV)
    call krome_set_user_crate(crate_0)
    call krome_set_user_Av(Av)

    dt = seconds_per_year * 1d4 !0.1 Myr initial time step
    t_tot = dt
    converged = .false.


    do i=1, rstep
      
      if (i .eq. 1) then
          !free-fall time, s
          tff = krome_get_free_fall_time(x(:))
          user_tff = tff         !store user tff. NEED THIS LINE FOR SOME WEIRD REASON
      endif

      ntot_val = sum(x(:))
      x(krome_idx_e) = krome_get_electrons(x(:))
      !Now do a renormalization of the number densities since the above changes number densities
      do j=1,krome_nmols
        x(j) = x(j) * ntot_val / sum(x(:))
      end do

      call krome_set_user_dissH2(dissH2)
      call krome_set_user_ionC(ionC)
      call krome_set_user_chiCO(chiCO)
      call krome_set_chiFUV(chiFUV)
      call krome_set_user_crate(crate_0)
      call krome_set_user_Av(Av)

      ni(krome_idx_Tgas) = Tgas

      !call equilibrium solver
      call krome_equilibrium_xT(x(:),Tgas,dt)
      Tdust = krome_get_Semenov_Tdust()

      !Reset Tgas to Tgas init
      Tgas = Tgasinit

      !avoid negative species
      do ii=1,krome_nmols
        n(ii) = max(x(ii),0d0)
      end do
      n(krome_idx_Tgas) = Tgas

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
                  "CONVERGED; zmetal = ", zmetal, " Tgas = ", Tgas, " t_tot/Myr = ", &
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
    rhogas = sum(x(:)*m(1:krome_nmols))
    write(22,'(99E17.8e3)') sum(x(:)),rhogas,Tgas,Tdust,zmetal,x(:)/sum(x(:))

    if (stop_next) exit

    !increase metallicity by a factor of 1.1 for the next bin
    zmetal = zmetal * 1.1
    !break when max density reached
    if (zmetal .gt. 1d0) then
      zmetal = 1d0
      stop_next = .true.
    endif
  end do

    !write(22, *)

  !Close files
  close(22)

  !say goodbye
  print *,"To plot in python:"
  print *,"ipython> run plot.py"
  print *,"That's all! have a nice day!"

end program test_krome_eqbm

