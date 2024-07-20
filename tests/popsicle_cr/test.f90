!################################################################
!This is a simple one-zone collapse test following
! the chemical and thermal evolution of a primordial cloud.
!The dynamics is described by the Larson-Penston-type
! similar solution and includes cooling and heating processes.
!For additional details look also to Omukai 2000 and
! the KROME paper.
!################################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_cooling
  implicit none
  integer,parameter::nz=8
  integer,parameter::rstep = 500000
  integer::i,unit,ios,jscale,jz,jz2
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),cools(krome_ncools)
  real*8::ntot,Tdust(krome_ndust),zs(nz)
  real*8::Av, NHtot, totheat, totcool, crate

  zs = (/0d0, 1d-6, 1d-5, 1d-4, 1d-3, 1d-2, 1d-1, 1d0/) !list of metallicities relative to solar

  !Cosmic ray ionization rate
  crate = 3d-17
  call krome_set_user_crate(crate)

  !output header
  write(22, '(A)', ADVANCE='NO') "#ntot Tgas Tdust"
  write(22, '(A)') trim(krome_get_names_header())

  write(31, '(A)', ADVANCE='NO') "#ntot Tgas sum(cools)"
  write(31, '(A)') trim(krome_get_cooling_names_header())

  !loop over size(zs)*2 so that every second loop is skipped, so that an empty line is created in the output fort.22 file
  !this line break in the output file can then be used to read in output for each zs separately
  do jz = 1,size(zs)*2

     jz2 = (jz+1)/2
     jscale = mod(jz,2)
     if(jscale==0) cycle
  

    !INITIAL CONDITIONS
    krome_redshift = 0d0    !redshift
    ntot = 0.1d0              !total density, cm-3
    Tgas = 3d2              !temperature, K

    call krome_set_zredshift(krome_redshift)
    call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
    call krome_set_metallicity(zs(jz2))

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
    x(KROME_idx_Cj)        = 0.927d-4*zs(jz2)*ntot !C is fully ionized
    x(KROME_idx_O)         = 3.568d-4*zs(jz2)*ntot !O is fully neutral

    call krome_init_dust_distribution(x(:),(1d0/162d0)*zs(jz2)) !scale the dust to gas ratio by the metallicity
    print *,krome_get_dust_distribution()
    call krome_set_Tdust((krome_redshift+1d0)*2.73d0)

    !set initial density
    dd = ntot

    !open file to write explore data
    open(newunit=unit,file="explore.dat",status="replace")

    print *,"solving..."
    print '(a5,3a11)',"step","n(cm-3)","Tgas(K)", "Tdust(K)"


    !print initial output
    Tdust = krome_get_Tdust()
    write(22,'(99E17.8e3)') dd,Tgas,Tdust(:),x(:)/dd

    !loop on density steps
    do i = 1,rstep

       !store old density
       dd1 = dd

       !free-fall time, s
       tff = krome_get_free_fall_time(x(:))
       user_tff = tff         !store user tff
       dtH = 0.01d0 * tff     !define time-step
       deldd = (dd/tff) * dtH !density increase
       dd = dd + deldd        !update density

       !rescale density
       x(:) = x(:)*dd/dd1

       !if you do not conserve electrons, the electron abundance will soon go to 0.00
       x(krome_idx_e) = krome_get_electrons(x(:))

       !set time-step
       dt = dtH

       !break when max density reached
       if(dd.gt.1d18) exit

       !important to scale the dust density as the gas density increases
       call krome_scale_dust_distribution(dd/dd1)

       !dust evaporation: dust is non existent at T > 1.5d3
       if(Tgas>1.5d3) call krome_scale_dust_distribution(0d0)

       Tdust = krome_get_Tdust()

       !dump cooling rates for Tgas going into the calculation
       n(1:krome_nmols) = x(:)
       n(KROME_idx_Tgas) = Tgas
       n(krome_nmols+krome_ndust+1:krome_nmols+2*krome_ndust) = Tdust
       cools(:) = get_cooling_array(n(:),Tgas)
       write(31,'(99E14.5e3)') dd, Tgas, sum(cools), cools(:)

       !solve the chemistry
       call krome(x(:),Tgas,dt)

       !print some output
       write(22,'(99E17.8e3)') dd,Tgas,Tdust(:),x(:)/dd
       if(mod(i,50)==0) then
          !totheat = krome_get_heating(x(:), Tgas)
          !totcool = krome_get_heating(x(:), Tgas)
          print '(I5,30E11.3)',i,dd,Tgas,Tdust(:)
          call krome_print_best_flux(x(:),Tgas,5)
          call krome_explore_flux(x(:),Tgas,unit,dd)
       end if
    end do
    write(22,*)
    write(31,*)
  end do

  !close explore data file
  close(unit)

  !say goodbye
  print *,"To plot in python:"
  print *,"ipython> run plot.py"
  print *,"That's all! have a nice day!"

end program test_krome
