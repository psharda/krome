! $Id: cooling.f90,v 1.40 2013/08/04 09:13:09 troels_h Exp $
!***********************************************************************
MODULE cooling_mod
  logical do_cool, do_radtrans
  !KROME: these variables are here for back-compatibilty
  real*8::T_MC
END MODULE cooling_mod

!***********************************************************************
SUBROUTINE read_cooling_namelist
  USE amr_commons, only: myid, chemistry
  USE cooling_mod
  implicit none
  namelist /cool/ do_cool,do_radtrans,chemistry
  rewind (1)
  read (1,cool)
  if (myid==1) write (*,cool)
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
  do_radtrans=.false.
  do_cool = .true.
  chemistry=.true.
  call read_cooling_namelist

  if(do_cool.or.chemistry) call krome_init()
  if (do_radtrans) call init_radiative_transfer

END SUBROUTINE init_cooling
