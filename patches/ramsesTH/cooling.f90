! $Id: cooling.f90,v 1.40 2013/08/04 09:13:09 troels_h Exp $
!***********************************************************************
MODULE cooling_mod
  logical :: do_cool, do_radtrans, do_oct_chemistry
  integer :: c_verbose ! verboseness level for cooling
  !KROME: these variables are here for back-compatibilty
  real*8::T_MC, Av_rho, crate, Tdust
END MODULE cooling_mod

!***********************************************************************
SUBROUTINE read_cooling_namelist
  USE amr_commons, only: myid, chemistry
  USE cooling_mod
  USE krome_user, only : krome_set_user_Tdust, krome_set_user_crate
  implicit none
  integer :: verbose  ! Local-var hack. A nice name in the namelist, but confilcts with global var "verbose"
  namelist /cool/ do_cool,do_radtrans,chemistry,verbose,crate,Av_rho,Tdust,do_oct_chemistry
  verbose = c_verbose ! Read module value
  rewind (1)
  read (1,cool)
  if (myid==1) write (*,cool)

  c_verbose = verbose ! write back to module value

  ! Make sure these are set. They are threadprivate.
  !$omp parallel
  call krome_set_user_Tdust(Tdust) ! Store it in the internal KROME structure
  call krome_set_user_crate(crate) ! Store it in the internal KROME structure
  !$omp end parallel

END SUBROUTINE read_cooling_namelist

!***********************************************************************
SUBROUTINE init_cooling
  USE amr_commons, only: print_id,chemistry,ncpu,myid
  USE cooling_mod
  USE krome_main
  USE krome_commons
  USE krome_user
  implicit none
  character(len=80):: id='$Id: cooling.f90,v 1.40 2013/08/04 09:13:09 troels_h Exp $'
!.......................................................................
  call print_id(id)
  if (ncpu>1) krome_mpi_rank=myid  ! Store the mpi rank in the corresponding krome variable
  do_radtrans = .false.
  do_cool     = .true.
  chemistry   = .true.
  do_oct_chemistry = .false.
  c_verbose   = 0
  crate       = 1.3e-17  ! Cosmic ray rate [s^-1]
  Tdust       = 10.      ! Temperature of the dust in Kelvin everywhere in the model (!)

  ! Normalise to Av = 1 for n ~ 1e3, and let it scale like 2/3 power.
  ! This is roughly correct according to Glover et al (astro-ph:1403.3530)
  Av_rho      = 0.001

  call read_cooling_namelist

  if(do_cool.or.chemistry) then
    !$omp parallel
    call krome_set_user_Av(1.0_8) ! Make sure AV is defined
    !$omp end parallel
    call krome_init()
  endif
  if (do_radtrans) call init_radiative_transfer

END SUBROUTINE init_cooling
