!##########################################################
!The following test provides a benchmark for
!gas-grain and grain-surface chemistry. 
!Initial conditions from Semenov+2010
!##########################################################
program test_krome

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_commons
  use krome_constants
  use krome_getphys

  integer,parameter::nx=krome_nmols
  integer::i
  real*8::x(nx),Tgas,t,dt,spy,xH,tout
  real*8::d2g,rho0,mmw,xdust,rhogas
  real*8::Av_opacity,CR_ion,adust,adust2
  real*8::mass(krome_nmols),dmass

  spy = 3600. * 24. * 365. !seconds per year
  Tgas = 1d1               !gas temperature (K)
  xH = 2d4                 !Hydrogen density xH  = n(H) + 2 n(H2)
  CR_ion = 1.3d-17         !cosmic rays ionisation rate in s^-1
  Av_opacity = 1d1         !exctintion coefficient
  d2g= 0.01                !dust to gas ratio
  rho0 = 3.0               !specific density for grain silicates
  mmw  = 1.43d0            !mean molecular weight
  adust = 1d-5             !0.1 micron average size for grains
  adust2=adust**2

  !user commons for opacity, CR rate, Tdust
  call krome_set_user_av(Av_opacity) !opacity Av (#)
  call krome_set_user_crflux(CR_ion)  !CR rate (1/s)
  call krome_set_user_Tdust(Tgas)         !set dust temperature (K)
  call krome_set_user_gsize(adust)     !set average dust size: 0.1 micron
  call krome_set_user_gsize2(adust2)     !set average dust size: 0.1 micron

  !initialise krome
  call krome_init()


  x(:) = 1.d-40
  !initial densities (Semenov2010)
  x(KROME_idx_He)  = 9d-2   * xH
  x(KROME_idx_H2)  = 0.5d0   * xH
  x(KROME_idx_Cj)  = 1.2d-4  * xH
  x(KROME_idx_N)   = 7.6d-5  * xH
  x(KROME_idx_O)   = 2.56d-4 * xH
  x(KROME_idx_Sj)  = 8.0d-8  * xH
  x(KROME_idx_Sij) = 8.0d-9  * xH
  x(KROME_idx_Fej) = 3d-9   * xH
  x(KROME_idx_Naj) = 2d-9   * xH
  x(KROME_idx_Mgj) = 7.0d-9  * xH
  x(KROME_idx_Pj)  = 2.0d-10 * xH
  x(KROME_idx_Clj) = 1.0d-9  * xH

  !calculate elctrons (neutral cloud)
  x(KROME_idx_e) = krome_get_electrons(x(:))

  !get dust density and other dust quantities
  !call krome_get_dust_variables(x(:),d2g,rho0)
  mass = krome_get_mass()
  rhogas = krome_get_rho(x(:))
  dmass=4./3d0*pi*rho0*adust**3
  xdust=d2g*rhogas/(dmass*1.43d0) !From Semenov 2010
  call krome_set_user_xdust(xdust)     !set average dust size: 0.1 micron
  x(KROME_idx_GRAIN0) =  xdust

  dt = 1d0*spy !time-step (s)
  t = 0d0 !initial time (s)

  !output header
  write(66,'(a)') "#time "//trim(krome_get_names_header())

  do
     print '(a10,E11.3,a3)',"time:",t/spy,"yr"
     call krome(x(:),Tgas,dt) !call KROME
     t = t + dt !increase time
     dt = max(1d2,t/3d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,x(:)/2d4
     if(t>1d9*spy) exit !exit when overshoot 1d8 years
  end do

  print *,"Output dump in fort.66"
  print *,"Default gnuplot: load 'plot.gps'"

  print *,"That's all! have a nice day!"

end program test_krome
