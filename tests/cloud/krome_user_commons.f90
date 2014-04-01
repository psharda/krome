module krome_user_commons
  implicit none

   ! *************************************************************
   !  This file has been generated with:
   !  krome 14.03 "Beastie Boyle" on 2014-04-01 18:43:06.
   ! 
   ! KROME is a nice and friendly chemistry package for a wide range of 
   !  astrophysical simulations. Given a chemical network (in CSV format) 
   !  it automatically generates all the routines needed to solve the kinetic 
   !  of the system, modelled as system of coupled Ordinary Differential 
   !  Equations. 
   !  It provides different options which make it unique and very flexible. 
   !  Any suggestions and comments are welcomed. KROME is an open-source 
   !  package, GNU-licensed, and any improvements provided by 
   !  the users is well accepted. See disclaimer below and GNU License 
   !  in gpl-3.0.txt.
   ! 
   !  more details in http://kromepackage.org/
   !  also see https://bitbucket.org/krome/krome_stable
   ! 
   ! Written and developed by Tommaso Grassi
   ! tommasograssi@gmail.com,
   ! University of Rome "Sapienza".
   ! 
   ! Co-developer Stefano Bovino
     ! sbovino@astro.physik.uni-goettingen.de
   ! Institut fuer Astrophysik, Goettingen.
   ! 
   ! Others (alphabetically): F.A. Gianturco, J.Prieto,
   ! D.R.G. Schleicher, D. Seifried, E. Simoncini, E. Tognelli.
   ! 
   ! 
   ! KROME is provided "as it is", without any warranty. 
   !  The Authors assume no liability for any damages of any kind 
   !  (direct or indirect damages, contractual or non-contractual 
   !  damages, pecuniary or non-pecuniary damages), directly or 
   !  indirectly derived or arising from the correct or incorrect 
   !  usage of KROME, in any possible environment, or arising from 
   !  the impossibility to use, fully or partially, the software, 
   !  or any bug or malefunction.
   !  Such exclusion of liability expressly includes any damages 
   !  including the loss of data of any kind (including personal data)
   ! *************************************************************

  
  !user can add here the common variables needed
  !for rate calculation (e.g. optical depth, CR rate, 
  !pressure, density, ...)
  real*8::tau,zrate,pah_size,gas_dust_ratio,krome_redshift
  real*8::krome_J21,user_tff
  real*8::krome_invdvdz !inverse of |velocity gradient| 1/abs(dv/dz)
  real*8::myCoe(4429)

contains 

  !user can add here the functions he/she needs for
  !rate calculations (Kooij funcion provided as example)
  function kooij(kalpha,kbeta,kgamma,Tgas)
    real*8::kooij,kalpha,kbeta,kgamma,Tgas
    kooij = kalpha*(Tgas/3d2)**kbeta*exp(-kgamma/Tgas)
  end function kooij


end module krome_user_commons
