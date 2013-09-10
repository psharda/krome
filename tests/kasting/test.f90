program test_krome

  !***************************
  !This test simulates the evolution of a planet atmosphere
  !for different atmospheric layers
  !***************************

  use krome_main
  use krome_user
  use krome_user_commons
  use krome_constants
 
  !nspec is shown during krome.py execution
  integer,parameter::nspec=57
  integer::i,j,ilay,ios
  real*8::x(nspec),dt,t,rho,Tgas,h,p,datar(59),xdummy

  open(22,file="layers_data",status="old")
  open(23,file="init_specs",status="old")
  read(22,*) !skip header
  !loop over layers
  do
     read(22,*,iostat=ios) ilay, h, p, Tgas, rho
     read(23,*) datar(:)
     if(ios.ne.0) exit
     x(krome_idx_H2O) = datar(1)
     x(krome_idx_O_1D) = datar(2)
     x(krome_idx_OH) = datar(3)
     x(krome_idx_H) = datar(4)
     x(krome_idx_H2) = datar(5)
     x(krome_idx_O) = datar(6)
     x(krome_idx_O2) = datar(7)
     x(krome_idx_O3) = datar(8)
     x(krome_idx_HO2) = datar(9)
     x(krome_idx_H2O2) = datar(11)
     x(krome_idx_N2) = datar(12)
     x(krome_idx_O_3P) = datar(13)
     x(krome_idx_CO) = datar(15)
     x(krome_idx_CO2) = datar(16)
     x(krome_idx_HCO) = datar(17)
     x(krome_idx_H2CO) = datar(18)
     x(krome_idx_CH4) = datar(19)
     x(krome_idx_CH3OOH) = datar(20)
     x(krome_idx_H3CO) = datar(21)
     x(krome_idx_N2O) = datar(22)
     x(krome_idx_HNO2) = datar(23)
     x(krome_idx_NO) = datar(24)
     x(krome_idx_HNO3) = datar(25)
     x(krome_idx_NO2) = datar(26)
     x(krome_idx_N) = datar(27)
     x(krome_idx_CH3) = datar(28)
     xdummy = datar(29) ! CH3O
     x(krome_idx_CH2_1) = datar(30)
     x(krome_idx_CH2_3) = datar(31)
     x(krome_idx_CH3O2) = datar(32)
     x(krome_idx_NO3) = datar(33)
     x(krome_idx_HO2NO2) = datar(34)
     x(krome_idx_CH3Cl) = datar(35)
     x(krome_idx_Cl) = datar(36)
     x(krome_idx_ClO) = datar(37)
     x(krome_idx_HCl) = datar(38)
     x(krome_idx_CH2Cl) = datar(39)
     x(krome_idx_ClONO2) = datar(40)
     x(krome_idx_NOCl) = datar(41)
     x(krome_idx_ClONO) = datar(42)
     x(krome_idx_Cl2) = datar(43)
     x(krome_idx_ClO2) = datar(44)
     x(krome_idx_HOCl) = datar(45)
     x(krome_idx_Cl2O2) = datar(46)
     x(krome_idx_N2O5) = datar(47)
     x(krome_idx_S) = datar(48)
     x(krome_idx_SO) = datar(49)
     x(krome_idx_SO2) = datar(50)
     x(krome_idx_H2S) = datar(51)
     x(krome_idx_HS) = datar(52)
     x(krome_idx_HSO3) = datar(53)
     x(krome_idx_SO3) = datar(54)
     x(krome_idx_H2SO4) = datar(55)
     x(krome_idx_S2) = datar(56)
     x(krome_idx_SO2_3) = datar(57)
     x(krome_idx_SO2_1) = datar(58)
     x(krome_idx_HSO) = datar(59)

     dt = seconds_per_year * 1d1 !initial time-step (s)
     t = 0.d0 !initial time (s)
     print *,"solving layer",ilay
     !loop over time to provide intermediate steps
     do 
        dt = dt * 2d0 !increase time-step (s)
        t = t + dt !increse time (s)
        call krome(x(:),Tgas,dt) !call KROME
        write(66,'(I5,99E12.3e3)') ilay,h,t/seconds_per_year, x(:) !write to file
        if(t>seconds_per_year*1d5) exit !break loop when overshoot 1e5 years
     end do
     write(77,'(I5,99E12.3e3)') ilay,h,t/seconds_per_year, x(:) !write to file
     write(66,*) !blank line at end of each loop
  end do
  print *,"done!"

  !bye bye
  print *,"that's all! have a nice day!"
  close(22)
  close(23)

end program test_krome
