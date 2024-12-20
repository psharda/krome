!################################################################
!One-zone test evolving temperatures and abundances to equilibrium, for a given fixed density, metallicity and radiation field.
!For additional details, see Sec. 7.1 Kim+23 ApJS
!Author: Shyam Menon (CCA/Rutgers, 2024); Piyush Sharda (Leiden, 2024)
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
  integer,parameter::nz=7
  integer,parameter::rstep = 500000
  integer::i,ii,ios,jscale,jz,jz2, dens_bins, zint
  real*8::rhogas,m(krome_nspec)
  real*8::tff,ertol,eatol,max_time,t_tot
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust,zs(nz),kk(krome_nrea),kkk(krome_nspec)
  real*8::Av,heats(krome_nheats),crate,NH,NHj,NH2
  real*8::ionH,dissH2,ionC,dissCO,chiFUV,t_cool
  logical::stop_next, converged
  character(len=20) :: filename, zint_str
  

  zs = (/1d-6, 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 1d0/) !list of metallicities relative to solar
  !zs = (/1d0/)

  !set chiFUV for photoreactions
  chiFUV = 1d0


  ! Switches to decide when equilibrium has been reached
  ertol = 1d-8  ! relative min change in a species
  max_time=seconds_per_year*5d10 ! max time we will be integrating for

  !loop over size(zs)*2 so that every second loop is skipped, so that an empty line is created in the output fort.22 file
  !this line break in the output file can then be used to read in output for each zs separately
  do jz = 1,size(zs)*2

    jz2 = (jz+1)/2
    jscale = mod(jz,2)
    if(jscale==0) cycle

    !Deduce filename from metallicity
    zint = int(log10(zs(jz2)))
    if(zint >= 0) then
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

    !CR rate proportional to FUV flux
    crate = 2d-16 * chiFUV
    print *, 'Initial crate: ', crate
    call krome_set_user_crate(crate)

    !switch to tell when to stop the calculation
    stop_next = .false.

    !Reset ertol for each metallicity (since it is changed at high density below)
    ertol = 1d-8

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
      x(KROME_idx_D)         = 3d-5*ntot

      call krome_set_Semenov_Tdust((krome_redshift+1d0)*2.73d0)

      !No shielding; Av=0.0
      Av = 0.0
      call krome_set_user_Av(Av)

      !set H2 dissociation reaction rate coeff
      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      ni(:) = n(:)

      !No LyC radiation field in the ISRF
      ionH = 0.0
      call krome_set_user_ionH(ionH)
      !set H2 dissociation reaction rate coeff
      dissH2 = 5.60d-11*chiFUV
      call krome_set_user_dissH2(dissH2)
      ionC = 3.1d-10*krome_get_user_is_metal()*chiFUV
      dissCO = 2.592d-10*krome_get_user_is_metal()*chiFUV
      call krome_set_user_ionC(ionC)
      call krome_set_user_dissCO(dissCO)

      !Start with a small initial timestep (10^4 yrs; based on the shortest cooling times in dense gas)
      dt = seconds_per_year * 1d4
      t_tot = dt
      converged = .false.

      !Higher densities, lower tolerance for convergence
      if(ntot .gt. 1.e5) ertol = 1d-6

      !loop on density steps
      do i=1, rstep

        if (i .eq. 1) then
          !free-fall time, s
          tff = krome_get_free_fall_time(x(:))
          user_tff = tff         !store user tff. NEED THIS LINE FOR SOME WEIRD REASON
        endif

        !if you do not conserve electrons, the electron abundance will soon go to 0.00
        x(krome_idx_e) = krome_get_electrons(x(:))

        Av = 0.0
        call krome_set_user_Av(Av)

        !set H ionization reaction rate coeff
        ionH = 0.0
        call krome_set_user_ionH(ionH)
        !set H2 dissociation reaction rate coeff
        !dissH2 = 5.60d-11*exp(-3.74*Av)*krome_fshield(n,Tgas)*chiFUV
        dissH2 = 5.60d-11*chiFUV
        call krome_set_user_dissH2(dissH2)
        ionC = 3.1d-10*krome_get_user_is_metal()*chiFUV
        dissCO = 2.592d-10*krome_get_user_is_metal()*chiFUV
        call krome_set_user_ionC(ionC)
        call krome_set_user_dissCO(dissCO)

        Tdust = krome_get_Semenov_Tdust()

        ni(krome_idx_Tgas) = Tgas

        !solve the chemistry and temperature evolution
        call krome_equilibrium_xT(x(:),Tgas,dt)

        !avoid negative species
        do ii=1,krome_nmols
          n(ii) = max(x(ii),0d0)
        end do
        n(krome_idx_Tgas) = Tgas

        kkk(1:krome_nmols) = krome_conserve(n(1:krome_nmols),ni(1:krome_nmols))
        n(:) = kkk(:)
        n(krome_idx_Tgas) = Tgas
        
        !Convergence test on temperature; also explored with including relative change in abundances, but did not make difference
        converged = abs(n(krome_idx_Tgas) - ni(krome_idx_Tgas)) / ni(krome_idx_Tgas) .le. ertol &
                     .or. t_tot .gt. max_time

        !Compute cooling time; t_cool = P/Lambda = nk_BT/Lambda; where Lambda is in erg cm^-3 s^-1
        t_cool = (sum(x(:)) * boltzmann_erg * Tgas)/(cooling(n(:),Tgas))

        if(.not. converged) then
          !Restrict timestep to 10% of cooling time for stability, and enforce at most factor 3 change in dt
          dt = MIN(0.1*t_cool,dt*3.0)
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
        print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
        print *, 'Tgas :', n(krome_idx_Tgas)
      end if

      m = get_mass()
      rhogas = sum(x(:)*m(1:krome_nmols))
      write(22,'(99E17.8e3)') sum(x(:)),rhogas,Tgas,Tdust,x(:)/sum(x(:))

      if (stop_next) exit

      !increase density by 2x for the next bin
      ntot = ntot * 1.1
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

  !say goodbye
  print *,"To plot in python:"
  print *,"ipython> run plot.py"
  print *,"That's all! have a nice day!"

end program test_krome_eqbm

