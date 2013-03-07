program test_krome

  use krome_main
  use krome_subs
  use krome_user
  use krome_user_commons

  integer,parameter::nspec=13,rstep = 5000
  integer::i
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(nspec),Tgas,rhogas,dt

  call krome_init()

  rhogas = 10.d0

  x(krome_idx_H) = 0.758d0
  x(krome_idx_E) = 0.121d0
  x(krome_idx_Hj) = 0.744d-3
  x(krome_idx_HE) = 0.250d-20
  x(krome_idx_HEj) = 0.250d-20
  x(krome_idx_HEjj) = 0.599d-1
  x(krome_idx_Hk) = 0.998d-20
  x(krome_idx_H2) = 0.759d-5
  x(krome_idx_H2j) = 0.759d-20


  Tgas = 0.826d2

  call get_infos(x(:), Tgas)

  dd = rhogas

  print *,"solving..."

  do i = 1,rstep

     dd1=dd

     !***CALCULATE THE FREE FALL TIME***!
     tff=sqrt(3.0d0*3.14d0/(32.0d0*6.67e-8*1.22d0*dd*1.67e-24))
     !dn / dt = n / t_ff(n)
     dtH=0.05d0*tff        !TIME-STEP
     
     user_tff = tff
     
     deldd=(dd/tff)*dtH
     dd=dd+deldd        !UPDATE DENSITY

     x(KROME_idx_H)    = x(KROME_idx_H)*dd/dd1  
     x(KROME_idx_E)    = x(KROME_idx_E)*dd/dd1  
     x(KROME_idx_Hj)   = x(KROME_idx_Hj)*dd/dd1 
     x(KROME_idx_HE)   = x(KROME_idx_HE)*dd/dd1
     x(KROME_idx_HEj)  = x(KROME_idx_HEj)*dd/dd1
     x(KROME_idx_H2)   = x(KROME_idx_H2)*dd/dd1
     x(KROME_idx_H2j)  = x(KROME_idx_H2j)*dd/dd1
     x(KROME_idx_Hk)   = x(KROME_idx_Hk)*dd/dd1
     x(KROME_idx_HEjj) = x(KROME_idx_HEjj)*dd/dd1
     print '(I5,99E11.3)',i,dd,Tgas
     write(22,*) dd,Tgas
     !call plot_cooling(x)

     dt = dtH 

     call krome(x(:),rhogas,Tgas,dt)

  end do

  print *,"done!"

  call get_infos(x(:), Tgas)

  print *,"that's all! have a nice day!"

end program test_krome
