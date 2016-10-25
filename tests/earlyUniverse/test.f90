!**********************************************************
! This test is aimed to reproduce the evolution
! of the primordial gas based on the Galli&Palla 1998 paper
!**********************************************************
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_constants

  integer,parameter::rstep = 5000
  integer::i
  real*8::dtH
  real*8::dd,dd1,rhoc
  real*8::x(krome_nmols),Tgas,dt,Trad
  real*8::ntot,dzeta,redshift
  real*8::aexp,Tcmb0
  real*8::zinit,zend,u

  !INITIALIZE KROME PARAMETERS
  call krome_init()

  !INITIAL CONDITIONS
  Tcmb0 = 2.73d0 !K
  rhoc = 9.19d-6 !h2/cm3
  zinit = 1d4 !2d3 !initial redshift
  aexp = 1d0 + zinit !1 + z factor
  ntot = krome_OmegaB * rhoc !total co-moving number density in 1/cm3
  dd = ntot * aexp**3 !proper initial number density
  Trad = Tcmb0 * aexp !initial radiation temperature in K
  Tgas = Trad !initial gas temperature in K

  !INTEGRATION PARAMETERS
  zend = 1d1 !final redshift
  dzeta = (zend - zinit) / (rstep - 1d0) !redshift step

  !species initialization in 1/cm3
  x(:) = 1d-40

  !INIT: FULLY IONIZED GAS IN NUMBER DENSITY
  x(KROME_idx_Hj)     = 0.924d0 * dd   !H+
  x(KROME_idx_HEjj)   = 0.076d0 * dd   !He++
  x(KROME_idx_Dj)     = 4.3d-5 * dd    !D+
  x(KROME_idx_E)      = krome_get_electrons(x)

  call krome_get_info(x(:),Tgas)

  print *,"solving..."
  print '(a7,4a11)',"step","z","n(cm-3)","Tgas(K)","Trad(K)"

  !header
  write(22,'(a)') "#z density Tgas Trad "//krome_get_names_header()

  !loop over the hydro time-step
  do i = 1,rstep

     dd1 = dd
     redshift = zinit + dzeta * (i-1)

     !H and D are computed using the approach
     ! of Peebles 1968 book. see eqn.(6.95)
     !note also that for redshift higher than 2.243d3
     ! the reactions involving everything except He
     ! are set to zero, by using the modifier found
     ! in the reaction network react_GP98
     if(redshift>2.243d3)then
        u = 1d1**(20.99d0 - log10(krome_Omegab * krome_hubble**2 &
             * (1d0 + redshift)**1.5) &
             - 2.5050d4 / (1d0 + redshift))
        x(krome_idx_H) = 0.924*dd / (2d0+u)
        x(krome_idx_D) = 4.3d-5*dd / (2d0+u)
     endif

     !set the internal redshift for krome
     call krome_set_zredshift(redshift)

     aexp = 1d0 + redshift  !update 1 + z

     !convert redshift step in timestep
     ! eq. (5) Galli&Palla 1998
     dtH = -dzeta / (krome_Hubble0 * aexp**2 &
          * sqrt(krome_Omega0 * redshift + 1d0))

     !evaluate radiation temperature
     Trad = Tcmb0 * (aexp)

     !set internal Trad
     call krome_set_user_Trad(Trad)

     !evaluate the new density
     dd = ntot * (aexp)**3

     !species update
     x(:) = x(:) * dd / dd1

     dt = dtH

     !solve the chemistry
     call krome(x(:),Tgas,dt)

     write(22,'(99E17.8e3)') aexp,dd,Tgas,Trad,x(:)/dd
     if(mod(i,100)==0) then
        print '(I7,99E11.3)',i,aexp,dd,Tgas,Trad
     end if
  end do

  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
