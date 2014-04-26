! $Id: cooling.f90,v 1.40 2013/08/04 09:13:09 troels_h Exp $
!***********************************************************************
MODULE cooling_mod
  logical :: do_cool, do_radtrans
  integer :: c_verbose ! verboseness level for cooling
  !KROME: these variables are here for back-compatibilty
  real*8::T_MC, Av_rho, crate
END MODULE cooling_mod

!***********************************************************************
SUBROUTINE read_cooling_namelist
  USE amr_commons, only: myid, chemistry
  USE cooling_mod
  implicit none
  use krome_user_commons, only : krome_crate
  integer :: verbose  ! Local-var hack. A nice name in the namelist, but confilcts with global var "verbose"
  namelist /cool/ do_cool,do_radtrans,chemistry,verbose,crate,Av_rho
  verbose = c_verbose ! Read module value
  rewind (1)
  read (1,cool)
  if (myid==1) write (*,cool)
  krome_crate = crate
  c_verbose = verbose ! write back to module value
END SUBROUTINE read_cooling_namelist

!***********************************************************************
SUBROUTINE init_cooling
  USE amr_commons, only: print_id,chemistry
  USE cooling_mod
  use krome_main
  implicit none
  character(len=80):: id='$Id: cooling.f90,v 1.40 2013/08/04 09:13:09 troels_h Exp $'
!.......................................................................
  call print_id(id)
  do_radtrans = .false.
  do_cool     = .true.
  chemistry   = .true.
  c_verbose   = 0
  crate       = 1.3e-17  ! Cosmic ray rate [s^-1]

  ! Normalise to Av = 1 for n ~ 1e3, and let it scale like 2/3 power.
  ! This is roughly correct according to Glover et al (astro-ph:1403.3530)
  Av_rho      = 0.001

  call read_cooling_namelist

  if(do_cool.or.chemistry) call krome_init()
  if (do_radtrans) call init_radiative_transfer

END SUBROUTINE init_cooling
