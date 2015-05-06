
!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  krome 14.08.dev "Beastie Boyle" on 2015-04-23 09:52:12.
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
  !  tommasograssi@gmail.com,
  !  Starplan Center, Copenhagen.
  !  Niels Bohr Institute, Copenhagen.
  !
  ! Co-developer Stefano Bovino
  !  sbovino@astro.physik.uni-goettingen.de
  !  Institut fuer Astrophysik, Goettingen.
  !
  ! Others (alphabetically): D. Galli, F.A. Gianturco, T. Haugboelle,
  !  J.Prieto, D.R.G. Schleicher, D. Seifried, E. Simoncini,
  !  E. Tognelli
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
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=0
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
contains
  
  !*******************
  function cooling(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(:),Tgas,cooling
    
    cooling = sum(get_cooling_array(n(:),Tgas))
    
  end function cooling
  
  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    
    !returns cooling in erg/cm3/s
    cools(:) = 0.d0
    
    cools(idx_cool_H2) = cooling_H2(n(:), Tgas)
    
    cools(idx_cool_atomic) = cooling_Atomic(n(:), Tgas)
    
    cools(idx_cool_dust) = cooling_dust(n(:), Tgas)
    
    cools(idx_cool_cont) = cooling_continuum(n(:), Tgas)
    
    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)
    
    get_cooling_array(:) = cools(:)
    
  end function get_cooling_array
  
  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom
    
    cooling_custom = 0d0
    
  end function cooling_custom
  
  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    use krome_subs
    use krome_commons
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)
    
    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3
    
    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return
    
    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0
    
    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))
    
    kpla = 1d1**(a0*y + a1) !fit density only
    
  end function kpla
  
  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas
    
    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90
    
    coolingChem = 0.d0
    
  end function coolingChem
  
  !**********************************
  function cooling_Continuum(n,Tgas)
    !cooling from continuum for a thin gas (no opacity)
    !see Omukai+2000 for details
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    real*8::n(:),Tgas,cooling_Continuum,kgas,rhogas
    real*8::lj,tau,beta,m(nspec)
    
    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3
    kgas = kpla(n(:),Tgas) !planck opacity cm2/g (Omukai+2000)
    lj = get_jeans_length(n(:), Tgas) !cm
    tau = lj * kgas * rhogas + 1d-40 !opacity
    beta = min(1.d0,tau**(-2)) !beta escape (always <1.)
    cooling_Continuum = 4.d0 * stefboltz_erg * Tgas**4 &
        * kgas * rhogas * beta !erg/s/cm3
    
  end function cooling_Continuum
  
  !*******************************
  function cooling_dust(n,Tgas)
    !cooling from dust in erg/cm3/s
    use krome_constants
    use krome_commons
    use krome_subs
    use krome_dust
    implicit none
    real*8::cooling_dust,n(:),Tgas
    real*8::rhogas,ljeans,be,ntot,vgas
    real*8::m(nspec),intCMB,fact
    integer::i
    fact = 0.5d0
    cooling_dust = 0d0
    m(:) = get_mass()
    Tgas = n(idx_Tgas)
    vgas = sqrt(kvgas_erg*Tgas) !thermal speed of the gas
    ntot = sum(n(1:nmols))
    rhogas = sum(n(1:nmols)*m(1:nmols))
    ljeans = get_jeans_length_rho(n(:),Tgas,rhogas)
    be = besc(n(:),Tgas,ljeans,rhogas)
    
    do i=1,ndust
      intCMB = get_dust_intBB(i,phys_Tcmb)
      cooling_dust = cooling_dust + (get_dust_intBB(i,n(nmols+ndust+i)) &
          - intCMB) * be * xdust(i) * krome_dust_asize2(i)
    end do
    cooling_dust = 4d0*pi*cooling_dust !erg/s/cm3
    return
    
    cooling_dust = dust_cooling !erg/s/cm3
    
  end function cooling_dust
  
  !*****************************
  !escape opacity for H2 cooling.
  ! courtesy of Kazu Omukai (2014)
  ! Einstein's A coefficients for spontaneous emission
  ! calculated by Turner, Kirby-Docken, & Dalgarno 1977, ApJS, 35, 281
  ! and the excitation energies for the levels of Borysow,
  ! Frommhold & Moraldi (1989), ApJ, 336, 495.
  function H2opacity_omukai(Tgas, n)
    use krome_commons
    use krome_subs
    implicit none
    real*8::H2opacity_omukai,Tgas,ntot,lTgas,lntot,n(:)
    
    ntot = sum(n(1:nmols))
    lTgas = log10(Tgas)
    lntot = log10(num2col(ntot,n(:)))
    
    H2opacity_omukai = 1d1**(fit_anytab2D(arrH2esc_ntot(:), &
        arrH2esc_Tgas(:), arrH2esc(:,:), xmulH2esc, &
        ymulH2esc,lntot,lTgas))
    
  end function H2opacity_omukai
  
  !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
  !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
  !UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
  !EACH SINGLE FUNCTION IS IN erg/s
  !FINAL UNITS = erg/cm3/s
  !*******************************
  function cooling_H2(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::n(:),Tgas
    real*8::temp,logt3,logt,cool,cooling_H2,T3
    real*8::LDL,HDLR,HDLV,HDL,fact
    real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
    real*8::dump63,dump14
    integer::i
    character*16::names(nspec)
    temp = max(Tgas, 1d1)
    
    T3 = temp * 1.d-3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0.d0
    
    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3
    
    !dumping function to extend 6e3 and 1e4 limits
    dump63 = 1d0/ (1d0 + exp(min((temp-1d4)*8d-4,3d2)))
    dump14 = 1d0/ (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
    
    !//H2-H
    if(temp>1d1 .and. temp<=1d2) then
      cool = cool +1.d1**(-16.818342D0 +3.7383713D1*logt3 &
          +5.8145166D1*logt32 +4.8656103D1*logt33 &
          +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
      cool = cool +1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
          -1.1332860D1*logt32 -2.7850082D1*logt33 &
          -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
      !note here that the limit has been extended from 6e3 to 1e6
    elseif(temp>1.d3 .and. temp<=1.d6) then
      cool = cool +1d1**(-2.4311209D1 +4.6450521D0*logt3 &
          -3.7209846D0*logt32 +5.9369081D0*logt33 &
          -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H) &
          * dump63
      
    end if
    
    !//H2-Hp, extended from 1e4 to 1e6
    if(temp>1.d1 .and. temp<=1.d4)  then
      cool = cool + 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
          +0.015391166d0*logt32 -0.23619985d0*logt33 &
          -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj) &
          * dump14
    end if
    
    !//H2-H2, limit extended from 6e3 to 1e6
    if(temp>1.d2 .and. temp<=1.d6) then
      cool = cool + 1d1**(-2.3962112D1 +2.09433740D0*logt3 &
          -.77151436D0*logt32 +.43693353D0*logt33 &
          -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) &
          * dump63
    end if
    
    !//H2-e
    if(temp>1d2 .and. temp<=5d2) then
      cool =  cool +1d1**(-2.1928796d1 + 1.6815730d1*logt3 &
          +9.6743155d1*logt32 +3.4319180d2*logt33 &
          +7.3471651d2*logt34 +9.8367576d2*logt35 &
          +8.0181247d2*logt36 +3.6414446d2*logt37 &
          +7.0609154d1*logt38)*n(idx_e)
      !note: limit extended from 1e4 to 1e6
    elseif(temp>5d2 .and. temp<1d6)  then
      cool = cool + 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
          +.93310622D0*logt32 +4.0406627d0*logt33 &
          -4.7274036d0*logt34 -8.8077017d0*logt35 &
          +8.9167183*logt36 + 6.4380698*logt37 &
          -6.3701156*logt38)*n(idx_e) &
          * dump63
    end if
    
    !//H2-He,  limit extended from 1e4 to 1e6
    if(temp>1.d1 .and. temp<=1d6) then
      cool =  cool + 1d1**(-2.3689237d1 +2.1892372d0*logt3&
          -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
          +.19191375d0*logt35)*n(idx_He)  &
          * dump63
    end if
    
    !check error
    if(cool>1.d30) then
      print *,"ERROR!!! H2 cooling <0 OR cooling >1.d30 erg/s/cm3"
      print *,"cool (erg/s/cm3): ",cool
      names(:) = get_names()
      do i=1,size(n)
        print '(I3,a18,E11.3)',i,names(i),n(i)
      end do
      stop
    end if
    
    cool = max(cool,0.d0)
    
    !this will avoid a division by zero and useless calculations
    if(cool==0.d0) then
      cooling_H2 = 0.d0
      return
    end if
    
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
    !UPDATED TO THE DATA BY GLOVER 2015, MNRAS
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
         3.e-24*exp(-0.51/t3)) !erg/s
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
    HDL  = (HDLR + HDLV) !erg/s
    HDL = HDL * dump14  
    
    LDL = cool !erg/s
    fact = HDL/LDL
    
    cooling_H2 = HDL*n(idx_H2)/(1.d0+(fact))  &
        * H2opacity_omukai(Tgas, n(:)) !erg/cm3/s
    
  end function cooling_H2
  
  !Atomic COOLING  Cen ApJS, 78, 341, 1992
  !UNITS = erg/s/cm3
  !*******************************
  function cooling_Atomic(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::Tgas,cooling_atomic,n(:)
    real*8::temp,T5,cool
    
    temp = max(Tgas,10.d0) !K
    T5 = temp/1.d5 !K
    cool = 0.d0 !erg/cm3/s
    
    !COLLISIONAL IONIZATION: H, He, He+, He(2S)
    cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)
    cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)
    cool = cool+ 4.95d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-6.31515d5/temp)*n(idx_e)*n(idx_Hej)
    cool = cool+ 5.01d-27*temp**(-0.1687)/(1.d0+sqrt(T5))&
        *exp(-5.5338d4/temp)*n(idx_e)**2*n(idx_Hej)
    
    !RECOMBINATION: H+, He+,He2+
    cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)
    cool = cool+ 1.55d-26*temp**(0.3647)*n(idx_e)*n(idx_Hej)
    cool = cool+ 3.48d-26*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hejj)
    
    !DIELECTRONIC RECOMBINATION: He
    cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
        *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)
    
    !COLLISIONAL EXCITATION:
    !H(all n), He(n=2,3,4 triplets), He+(n=2)
    cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)
    cool = cool+ 9.1d-27*temp**(-.1687)/(1.d0+sqrt(T5))&
        *exp(-1.3179d4/temp)*n(idx_e)**2*n(idx_Hej)
    cool = cool+ 5.54d-17*temp**(-.397)/(1.d0+sqrt(T5))&
        *exp(-4.73638d5/temp)*n(idx_e)*n(idx_Hej)
    
    cooling_atomic = max(cool, 0.d0)  !erg/cm3/s
    
  end function cooling_Atomic
  
  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab
    
    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)
    
  end subroutine mylin2
  
  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)
    
    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if
    
    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if
    
    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)
    
  end subroutine mylin3
  
  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(33,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_H2 = cooling_H2(n(:),Tgas)
      cool_atomic = cooling_atomic(n(:),Tgas)
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(33)
    print *,"done!"
    
  end subroutine plot_cool
  
  !***********************************
  !routine to dump cooling in unit nfile
  subroutine dump_cool(n,Tgas,nfile)
    implicit none
    real*8::Tgas,n(:),cool(9)
    integer::nfile
    
    cool(:) = 0.d0
    
    cool(1) = cooling_H2(n(:), Tgas)
    
    cool(3) = cooling_atomic(n(:), Tgas)
    
    cool(7) = cooling_dust(n(:), Tgas)
    
    write(nfile,'(99E14.5e3)') Tgas,sum(cool),cool(:)
    
  end subroutine dump_cool
  
end module KROME_cooling
