program test_krome

  use krome_main
  use krome_subs
  use krome_user
  use krome_user_commons
  use krome_constants

  integer,parameter::nx=452
  real*8::x(nx),Tgas,t,dt,spy,xH

  spy = seconds_per_year
  Tgas = 1d1 !gas temperature (K)
  xH = 2d4 !Hydrogen density

  !user commons for opacity and CR rate
  tau = 1d1 !opacity Av (#)
  zrate = 1.3d-17 !CR rate (1/s)
  gas_dust_ratio = 7.57d11 !gas/dust
  pah_size = 4d-8 !cm

  call krome_init()

  x(:) = 1.d-20
  !initial densities (model EA2 Wakelam+Herbst 2008)
  x(KROME_idx_H2)  = 0.5d0   * xH 
  x(KROME_idx_He)  = 9.d-2   * xH 
  x(KROME_idx_N)   = 7.6d-5  * xH 
  x(KROME_idx_O)   = 2.56d-4 * xH 
  x(KROME_idx_Cj)  = 1.2d-4  * xH 
  x(KROME_idx_Sj)  = 1.5d-5  * xH 
  x(KROME_idx_Sij) = 1.7d-6  * xH
  x(KROME_idx_Fej) = 2.d-7   * xH 
  x(KROME_idx_Naj) = 2.d-7   * xH 
  x(KROME_idx_Mgj) = 2.4d-6  * xH 
  x(KROME_idx_Clj) = 1.8d-7  * xH 
  x(KROME_idx_Pj)  = 1.17d-7 * xH
  x(KROME_idx_Fj)  = 1.8d-8  * xH

  !calculate elctrons (neutral cloud)
  x(KROME_idx_e) = x(KROME_idx_Cj) + x(KROME_idx_Sj) + x(KROME_idx_Sij) &
       + x(KROME_idx_Fej) + x(KROME_idx_Naj) + x(KROME_idx_Mgj) &
       + x(KROME_idx_Clj) + x(KROME_idx_Pj) + x(KROME_idx_Fj)   

  dt = 1d2*spy !time-step (s)
  t = 0.d0 !initial time (s)
  write(66,'(999E12.3e3)') t/spy,x(:)
  do
     print '(a10,E11.3,a3)',"time:",t/spy,"yr"
     call krome(x(:),Tgas,dt) !call KROME
     t = t + dt !increase time
     dt = max(1d2,t/3.d0) !increase time-step
     write(66,'(999E12.3e3)') t/spy,x(:)
     if(t>1d8*spy) exit !exit when overshoot 1d8 years
  end do
 
  print *,"Output dump in fort.66"
  print *,"Default gnuplot: load 'plot.gps'"
  
  print *,"That's all! have a nice day!"

end program test_krome
