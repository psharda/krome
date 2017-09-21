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
  implicit none
  integer,parameter::rstep = 500000
  integer::i,unit,ios
  real*8::dtH,deldd
  real*8::tff,dd,dd1
  real*8::x(krome_nmols),Tgas,dt
  real*8::ntot

  !INITIAL CONDITIONS
  krome_redshift = 15d0    !redshift
  ntot = 1d0               !total density, cm-3
  Tgas = 1d2               !temperature, K

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

  !set initial density
  dd = ntot

  !open file to write explore data
  open(newunit=unit,file="explore.dat",status="replace")

  print *,"solving..."
  print '(a5,2a11)',"step","n(cm-3)","Tgas(K)"

  !output header
  write(22,*) "#ntot Tgas "//trim(krome_get_names_header())

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

     !set time-step
     dt = dtH

     !break when max density reached
     if(dd.gt.1d18) exit

     !solve the chemistry
     call krome(x(:),Tgas,dt)

     !print some output
     write(22,'(99E17.8e3)') dd,Tgas,x(:)/dd
     if(mod(i,100)==0) then
        print '(I5,99E11.3)',i,dd,Tgas
        call krome_print_best_flux(x(:),Tgas,5)
        call krome_explore_flux(x(:),Tgas,unit,dd)
     end if

  end do

  !close explore data file
  close(unit)

  !say goodbye
  print *,"To plot type in gnuplot:"
  print *,"gnuplot> load 'plot.gps'"
  print *,"That's all! have a nice day!"

end program test_krome
