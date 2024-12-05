!################################################################
!One-zone test where we fix temperature, solve the abundances to equilibrium, and obtain the cooling rates for this state
!This is a simple test to check the cooling rates for a given temperature; target figure would be something like Fig 3. of Kim+23
!Author: Shyam Menon (CCA/Rutgers, 2024)
!Email: smenon@flatironinstitute.org
!################################################################

program test_cooling
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
    integer::i,ii,ios,jscale,jz2, dens_bins, zint
    real*8::rhogas,m(krome_nspec)
    real*8::tff,ertol,eatol,max_time,t_tot
    real*8::x(krome_nmols),Tgas,dt,n(krome_nspec),ni(krome_nspec),cools(krome_ncools)
    real*8::ntot,Tdust,zs(nz),kk(krome_nrea),kkk(krome_nspec)
    real*8::Av,heats(krome_nheats),crate,NH,NHj,NH2
    real*8::ionH,dissH2,ionC,dissCO,chiFUV,t_cool,logTgas, f1, f2, Tt1,Tt2
    logical::stop_next, converged
    character(len=20) :: filename, zint_str

    zs = (/1d0/)
    jz2 = 1 !Only the first element

    !INITIAL CONDITIONS
    krome_redshift = 0d0    !redshift

    call krome_set_zredshift(krome_redshift)
    call krome_set_Tcmb(2.73d0*(krome_redshift+1d0))
    call krome_set_metallicity(zs(jz2))
    call krome_set_user_is_metal(1d0)

    !initialize KROME (mandatory)
    call krome_init()

    !CR rate proportional to FUV flux
    chiFUV = 1.0
    crate = 2d-16 * chiFUV
    print *, 'Initial crate: ', crate
    call krome_set_user_crate(crate)
    call krome_set_chiFUV(chiFUV)

    
    !Density is fixed
    ntot = 1.0
    !Initial temperature
    Tgas = 4936

    x(:) = 1d-40 !Default

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
    Av = 0.0
    call krome_set_user_Av(Av)

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

    tff = krome_get_free_fall_time(x(:))
    user_tff = tff

    open(unit=45,file="coolingfunc.dat",status='replace',action='write')

    do while(Tgas<1.e8)

      n(1:krome_nmols) = x(:)
      n(KROME_idx_Tgas) = Tgas
      ni(:) = n(:)

      x(krome_idx_e) = krome_get_electrons(x(:))

      !Solve abundances to equilibrium
      call krome_equilibrium(x(:),Tgas,1)

      x(krome_idx_e) = krome_get_electrons(x(:))


      !avoid negative species
      do ii=1,krome_nmols
        n(ii) = max(x(ii),0d0)
      end do
      n(krome_idx_Tgas) = Tgas

      !Compute smoothing factors for obtaining individual cooling rates
      !f2 becomes S in equation 40 of Kim+2023 ApJS, 264, 10
      Tt1 = 2d4
      Tt2 = 3.5d4
      f2 = 1d0/(1d0 + exp(-1d1*(Tgas - 0.5d0*(Tt1 + Tt2))/ (Tt2 - Tt1)))
      !f1 becomes (1-S) in equation 40 of Kim+2023 ApJS, 264, 10
      f1 = 1d0 - f2

      print * , "Tgas, Total Cooling: ", Tgas, cooling(n,Tgas)/ntot**2, "erg cm^3 s^-1"
      write(45,'(8E14.5e3)') ntot, Tgas, cooling(n,Tgas)/ntot**2, f1 * cooling_Nebular(n, Tgas)/ntot**2, f2 * cooling_Z_CIEGF(n, Tgas)/ntot**2,&
                                                            cooling_Atomic(n(:), Tgas)/ntot**2, f1 * cool_DustGRREC(n,Tgas)/ntot**2, f1 * ( cooling_Z(n(:), Tgas) - cooling_Z(n(:),2.73d0*(krome_redshift+1d0)) )/ntot**2

      !Make sure output is written before moving on
      flush(45)
      
      Tgas = Tgas*1.1

    end do

    print *, "Done"

end program test_cooling
