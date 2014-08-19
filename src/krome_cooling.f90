module KROME_cooling
#KROME_header
  integer,parameter::coolTab_n=int(1e2)
#KROME_nZrate
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
#KROME_escape_vars
#KROME_coolingZ_popvars
contains

  !*******************
  function cooling(n,Tgas)
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

#IFKROME_useCoolingH2
    cools(idx_cool_H2) = cooling_H2(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingH2GP
    cools(idx_cool_H2GP) = cooling_H2GP(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingAtomic
    cools(idx_cool_atomic) = cooling_Atomic(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingHD
    cools(idx_cool_HD) = cooling_HD(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingZ
    cools(idx_cool_Z) = cooling_Z(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingdH
    cools(idx_cool_dH) = cooling_dH(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDust
    cools(idx_cool_dust) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCompton
    cools(idx_cool_compton) = cooling_compton(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCIE
    cools(idx_cool_CIE) = cooling_CIE(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingContinuum
    cools(idx_cool_cont) = cooling_continuum(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingExpansion
    cools(idx_cool_exp) = cooling_expansion(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingFF
    cools(idx_cool_ff) = cooling_ff(n(:), Tgas)
#ENDIFKROME

    get_cooling_array(:) = cools(:)

    !remove the comment below to write cooling contributions to fort.44
    !write(44,'(99E17.8e3)') sum(n(1:nmols)),Tgas,cools(:)

    !gnuplot command (every n=100, and m=1 for density or m=2 for 
    ! temperature following the write(44,...) command above) 
    !plot 'fort.44' u m:3 every n w l t "H2",\
    ! '' u m:4 every n w l t "H2GP",\
    ! '' u m:5 every n w l t "atomic",\
    ! '' u m:6 every n w l t "HD",\
    ! '' u m:7 every n w l t "Z",\
    ! '' u m:8 every n w l t "enthalpy",\
    ! '' u m:9 every n w l t "dust",\
    ! '' u m:10 every n w l t "compton",\
    ! '' u m:11 every n w l t "CIE",\
    ! '' u m:12 every n w l t "Cont",\
    ! '' u m:13 every n w l t "Exp"

  end function get_cooling_array


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
    real*8::a0,a1

    rhogas = sum(n(:)*get_mass()) !g/cm3

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


!!$  !******************************
!!$  function cooling_FB(n,Tgas)
!!$    !free-bound cooling, as in Mewe+1986, sect2.2
!!$    ! using the implementation as in Chianti 
!!$    ! procedure fb_rad_loss.pro (erg/cm3/s)
!!$    use krome_commons
!!$    implicit none
!!$    real*8::coolingFB,n(:),Tgas,T6,cool,gfb,ln1
!!$    real*8::zetas()
!!$
!!$
!!$    T6 = Tgas*1d-6
!!$    cool = 0d0
!!$
!!$#KROME_FB_cooling_data
!!$
!!$
!!$    n0 = 
!!$    z0 = sqrt(ionpot/eH)*n0
!!$
!!$    e0 = ionpot/8.06554d3/1d3
!!$    en1 = eH*ion**2/(n0+1d0**2)/8.06554d3/1d3
!!$
!!$    l0 = 1d8/(ionpot-ionpot0)
!!$    ln1 = 1d8/eH/ion**2*(n0+1d0)**2
!!$    zeta0 = zetas(iz-ion)
!!$
!!$    precool = 2.051d-22*sqrt(T6)/143.9d0
!!$    pregf = 0.1578/T6*nion*ieq
!!$
!!$    f2 = 0.9*zeta0*z0**4/n0**5*exp(0.1578*z0**2/n0**2*T6)
!!$    gf = pregf*f2 
!!$    cool = cool + precool * exp(-143.9/l0/T6)
!!$
!!$    f2 = 0.42d0*n0**(-1.5)*ion**4 * exp(.1578*(ion/(n0+1d0))**2*T6)
!!$    gfb = pregf*f2
!!$    cool = cool + precool * gfb * exp(-143.9/ln1/T6)
!!$    
!!$    cooling_FB = cool
!!$
!!$  end function cooling_FB

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

#IFKROME_useCoolingContinuum
  !**********************************
  function cooling_Continuum(n,Tgas)
    !cooling from continuum for a thin gas (no opacity)
    !see Omukai+2000 for details
    use krome_commons
    use krome_constants
    use krome_subs
    implicit none
    real*8::n(:),Tgas,cooling_Continuum,kgas,rhogas
    real*8::lj,tau,beta
    rhogas = sum(n(:)*get_mass()) !g/cm3
    kgas = kpla(n(:),Tgas) !planck opacity cm2/g (Omukai+2000)
    lj = get_jeans_length(n(:), Tgas) !cm
    tau = lj * kgas * rhogas + 1d-40 !opacity
    beta = min(1.d0,tau**(-2)) !beta escape (always <1.)
    cooling_Continuum = 4.d0 * stefboltz_erg * Tgas**4 &
         * kgas * rhogas * beta !erg/s/cm3

  end function cooling_Continuum
#ENDIFKROME


#IFKROME_useCoolingCIE
  !*******************************
  function cooling_CIE(n, inTgas)
    !CIE cooling: fit from Ripamponti&Abel2004 (RA04) data
    ! The fit is valid from 100K-1e6K.
    ! Original data are from 400K to 7000K.
    ! We extrapolated data under 400K and fitted from 100K to 10**2.95 K.
    ! Data from 10**2.95 K to 1e5K are fitted analogously.
    ! Above 1e5 we employ a cubic extrapolation.
    use krome_commons
    use krome_constants
    real*8::cooling_CIE,n(:),Tgas,inTgas
    real*8::x,x2,x3,x4,x5
    real*8::a0,a1,a2,a3,a4,a5
    real*8::b0,b1,b2,b3,b4,b5
    real*8::cool,tauCIE,logcool

    !under 1e-12 1/cm3 cooling is zero
    if(n(idx_H2)<1d-12) then
       cooling_CIE = 0.d0
       return
    end if

    !temperature limit
    Tgas = max(inTgas, phys_Tcmb) 

    !prepares variables
    x = log10(Tgas)
    x2 = x*x
    x3 = x2*x
    x4 = x3*x
    x5 = x4*x

    cool = 0.d0
    !outside boundaries below cooling is zero
    logcool = -1d99

    !evaluates fitting functions
    if(x>2.d0 .and. x<2.95d0) then
       a0 = -30.3314216559651d0
       a1 = 19.0004016698518d0
       a2 = -17.1507937874082d0 
       a3 = 9.49499574218739d0 
       a4 = -2.54768404538229d0 
       a5 = 0.265382965410969d0
       logcool = a0 + a1*x + a2*x2 + a3*x3 +a4*x4 +a5*x5 
    elseif(x.GE.2.95d0 .and. x<5.d0) then
       b0 = -180.992524120965d0 
       b1 = 168.471004362887d0 
       b2 = -67.499549702687d0 
       b3 = 13.5075841245848d0 
       b4 = -1.31983368963974d0 
       b5 = 0.0500087685129987d0
       logcool = b0 + b1*x + b2*x2 + b3*x3 +b4*x4 +b5*x5
    elseif(x.GE.5.d0) then
       logcool = 3.d0 * x - 21.2968837223113 !cubic extrapolation
    end if

    !opacity according to RA04
    tauCIE = (n(idx_H2) * 1.4285714e-16)**2.8 !note: 1/7e15 = 1.4285714e-16
    cool = p_mass * 1d1**logcool !erg*cm3/s

    cooling_CIE = cool * min(1.d0, (1.d0-exp(-tauCIE))/tauCIE) &
         * n(idx_H2) * sum(n(1:nmols)) !erg/cm3/s

  end function cooling_CIE
#ENDIFKROME


#IFKROME_useCoolingExpansion
  !*******************************
  function cooling_expansion(n, Tgas)
    !R'/R expansion cooling erg/cm3/s from Galli&Palla 1998
    use krome_user_commons
    use krome_commons
    use krome_constants
    real*8::cooling_expansion,n(:),Tgas
    real*8::ntot

    ntot=sum(n(1:nmols))
    !note that user_* must be provided by the user 
    ! in the reaction file using @common: and initialized
    ! using the interface subroutine
    cooling_expansion = 3.d0*ntot*boltzmann_erg*Tgas*hubble0 &
         * (1.d0 + phys_zredshift) &
         * sqrt(omega0 * phys_zredshift + 1.d0) !erg/s/cm3

  end function cooling_expansion
#ENDIFKROME

#IFKROME_useCoolingCompton
  !*******************************
  function cooling_compton(n, Tgas)
    !compton cooling erg/cm3/s from Cen1992
    use krome_user_commons
    use krome_commons
    real*8::cooling_compton,n(:),Tgas

    !note that redhsift is defined in krome_user_commons and 
    ! must be provided by the user
    cooling_compton = 5.65d-36 * (1.d0 + phys_zredshift)**4 &
         * (Tgas - 2.73d0 * (1.d0 + phys_zredshift)) * n(idx_e) !erg/s/cm3

  end function cooling_compton
#ENDIFKROME

#IFKROME_useCoolingDust
  !*******************************
  function cooling_dust(n,Tgas)
    !cooling from dust in erg/cm3/s
    use krome_dust
    use krome_commons
    use krome_constants
    implicit none
    real*8::cooling_dust,n(:),Tgas,cool,ntot
    integer::i,idust

    !total gas density 1/cm3
    ntot = sum(n(1:nmols))

    cool = 0.d0 !cooling in erg/s/cm3
    !loop on dust to evaluate cooling following Hollenbach and McKee 1979
    do i=nmols+1,nmols+ndust
       idust = i - nmols !index of dust
       cool = cool + dustCool(krome_dust_asize2(idust), n(i), Tgas,&
            krome_dust_T(idust), ntot)
    end do

    cooling_dust = cool !erg/s/cm3

  end function cooling_dust
#ENDIFKROME


#IFKROME_useCoolingdH
  !*******************************
  function cooling_dH(n,Tgas)
    !cooling from reaction enthalpy erg/s/cm3
    use krome_commons
    implicit none
    real*8::cooling_dH,cool,n(:),Tgas,small,nmax
    real*8::logT,lnT,Te,lnTe,T32,t3,invT,invTe,sqrTgas,invsqrT32,sqrT32
#KROME_vars

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = #KROME_small

    logT = log10(Tgas) !log10 of Tgas (#)
    lnT = log(Tgas) !ln of Tgas (#)
    Te = Tgas*8.617343d-5 !Tgas in eV (eV)
    lnTe = log(Te) !ln of Te (#)
    T32 = Tgas/3.d2 !Tgas/(300 K) (#)
    t3 = T32 !alias for T32 (#)
    invT = 1.d0/Tgas !inverse of T (1/K)
    invTe = 1.d0/Te !inverse of T (1/eV)
    sqrTgas = sqrt(Tgas) !Tgas rootsquare (K**0.5)
    invsqrT32 = 1.d0/sqrt(T32)
    sqrT32 = sqrt(T32)

    cool = 0.d0

#KROME_rates
#KROME_dH_cooling

    cooling_dH = cool    

  end function cooling_dH
#ENDIFKROME


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

#IFKROME_useCoolingH2GP
  !*******************************
  function cooling_H2GP(n, Tgas)
    !cooling from Galli&Palla98
    use krome_commons
    use krome_subs
    real*8::n(:),Tgas, tm, logT
    real*8::cooling_H2GP,T3
    real*8::LDL,HDLR,HDLV,HDL,fact

    tm = max(Tgas, 13.0d0)    ! no cooling below 13 Kelvin
    tm = min(Tgas, 1.d5)      ! fixes numerics
    logT = log10(tm)
    T3 = tm * 1.d-3

    !low density limit in erg/s
    LDL = 1.d1**(-103.d0+97.59d0*logT-48.05d0*logT**2&
         +10.8d0*logT**3-0.9032d0*logT**4)*n(idx_H)

    !this will avoid a division by zero and useless calculations
    if(LDL==0d0) then
       cooling_H2GP = 0.d0
       return
    end if

    !high density limit
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
         3.e-24*exp(-0.51/t3)) !erg/s
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
    HDL  = HDLR + HDLV !erg/s

    !TO AVOID DIVISION BY ZERO
    fact = HDL/LDL !dimensionless

    cooling_H2GP = HDL*n(idx_H2)/(1.d0+(fact)) #KROME_H2opacity !erg/cm3/s


  end function cooling_H2GP
#ENDIFKROME

#IFKROME_useCoolingH2
    
  !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
  !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
  !EACH SINGLE FUNCTION IS IN erg/s
  !FINAL UNITS = erg/cm3/s
  !*******************************
  function cooling_H2(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::n(:),Tgas
    real*8::temp,logt3,logt,cool,cooling_H2,T3
    real*8::LDL,HDLR,HDLV,HDL,fact
    real*8::logt32,logt33,logt34,logt35,dump63,dump14
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
       cool = cool + 1d1**(-2.1716699D1 +1.3865783D0*logt3 &
            -0.37915285D0*logt32 +0.11453688D0*logt33 &
            -0.23214154D0*logt34 +0.058538864D0*logt35)*n(idx_Hj) &
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
    if(temp>1d1 .and. temp<=2d2) then
       cool =  cool +1d1**(-3.4286155D1 -4.8537163D1*logt3 &
            -7.7121176D1*logt32 -5.1352459D1*logt33 &
            -1.5169150D1*logt34 -.98120322D0*logt35)*n(idx_e)
       !note: limit extended from 1e4 to 1e6
    elseif(temp>2d2 .and. temp<1d6)  then
       cool = cool + 1d1**(-2.2190316D1 +1.5728955D0*logt3 &
            -.213351D0*logt32 +.96149759D0*logt33 &
            -.91023195D0*logt34 +.13749749D0*logt35)*n(idx_e) &
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

    !high density limit from HM79, GP98 
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL 
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
         3.e-24*exp(-0.51/t3))
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3))
    HDL  = HDLR + HDLV

    LDL = cool !erg/s
    fact = HDL/LDL 

    cooling_H2 = HDL*n(idx_H2)/(1.d0+(fact))  #KROME_H2opacity !erg/cm3/s

  end function cooling_H2
#ENDIFKROME


#IFKROME_useCoolingAtomic
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
#ENDIFKROME

#IFKROME_useCoolingFF

  !**************************
  !free-free cooling (bremsstrahlung for all ions)
  ! using mean Gaunt factor value (Cen+1992)
  function cooling_ff(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(:),Tgas,cool,cooling_ff,gaunt_factor,bms_ions

    gaunt_factor = 1.5d0 !mean value

    !BREMSSTRAHLUNG: all ions
#KROME_brem_ions
    cool = 1.42d-27*gaunt_factor*sqrt(Tgas)&
         *bms_ions*n(idx_e)

    cooling_ff = max(cool, 0.d0)  !erg/cm3/s

  end function cooling_ff
#ENDIFKROME

#IFKROME_useCoolingHD
  !HD COOLING LIPOVKA ET AL. MNRAS, 361, 850, (2005)
  !UNITS=erg/cm3/s
  !*******************************
  function cooling_HD(n, inTgas)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer, parameter::ns=4
    real*8::cooling_HD
    real*8::n(:),Tgas,logTgas,lognH,inTgas
    real*8::c(0:ns,0:ns),logW,W,dd,lhj

    !default HD cooling value
    cooling_HD = 0.0d0 !erg/cm3/s

    !this function does not have limits on density
    ! and temperature, even if the original paper do.
    ! However, we extrapolate the limits.

    !exit on low temperature
    if(inTgas<phys_Tcmb) return
    !extrapolate higher temperature limit
    Tgas = min(inTgas,1d4)

    !calculate density
    dd = n(idx_H) !sum(n(1:nmols))
    !exit if density is out of Lipovka bounds (uncomment if needed)
    !if(dd<1d0 .or. dd>1d8) return

    !extrapolate density limits
    dd = min(max(dd,1d-2),1d10)

    !POLYNOMIAL COEFFICIENT: TABLE 1 LIPOVKA 
    c(0,:) = (/-42.56788d0, 0.92433d0, 0.54962d0, -0.07676d0, 0.00275d0/)
    c(1,:) = (/21.93385d0, 0.77952d0, -1.06447d0, 0.11864d0, -0.00366d0/)
    c(2,:) = (/-10.19097d0, -0.54263d0, 0.62343d0, -0.07366d0, 0.002514d0/)
    c(3,:) = (/2.19906d0, 0.11711d0, -0.13768d0, 0.01759d0, -0.00066631d0/)
    c(4,:) = (/-0.17334d0, -0.00835d0, 0.0106d0, -0.001482d0, 0.00006192d0/)

    logTgas = log10(Tgas)
    lognH   = log10(dd)

    !loop to compute coefficients
    logW = 0.d0
    do j = 0, ns
       lHj = lognH**j
       do i = 0, ns
          logW = logW + c(i,j)*logTgas**i*lHj !erg/s
       enddo
    enddo

    W = 10.d0**(logW)
    cooling_HD = W * n(idx_HD) !erg/cm3/s


  end function cooling_HD
#ENDIFKROME

#IFKROME_useCoolingZ
  !*********************************************
  !function for linear interpolation of f(x), using xval(:)
  ! and the corresponding yval(:) as reference values
  ! note: slow function, use only for initializations
  function flin(xval,yval,x)
    implicit none
    real*8::xval(:),yval(:),x,flin
    integer::i,n
    logical::found
    found = .false.
    n = size(xval)
    x = max(x,xval(1)) !set lower bound
    x = min(x,xval(n)) !set upper bound
    !loop to find interval (slow)
    do i=2,n
       if(x.le.xval(i)) then
          !linear fit
          flin = (yval(i) - yval(i-1)) / (xval(i) - xval(i-1)) * &
               (x - xval(i-1)) + yval(i-1)
          found = .true. !found flag
          exit
       end if
    end do
    if(.not.found) flin = yval(n)

  end function flin

  !************************
  !dump the level populations in a file
  subroutine dump_cooling_pop(Tgas,nfile)
    implicit none
    integer::nfile,i
    real*8::Tgas
    
#KROME_popvar_dump
    write(nfile,*)

  end subroutine dump_cooling_pop

  !***********************
  !metal cooling as in Maio et al. 2007
  ! loaded from data file 
  function cooling_Z(n,inTgas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:), inTgas, cool, cooling_Z, k(nZrate), Tgas

    Tgas = inTgas
    k(:) = coolingZ_rate_tabs(Tgas)

    cool = 0d0
#KROME_coolingZ_call_functions

    cooling_Z = cool * boltzmann_erg

  end function cooling_Z

  !********************************
  function coolingZ_rates(inTgas)
    use krome_commons
    use krome_subs
    implicit none
    real*8::inTgas,coolingZ_rates(nZrate),k(nZrate)
    real*8::Tgas,invT,logTgas
    integer::i
#KROME_coolingZ_declare_custom_vars

    Tgas = inTgas
    invT = 1d0/Tgas
    logTgas = log10(Tgas)

#KROME_coolingZ_custom_vars

#KROME_coolingZ_rates

    coolingZ_rates(:) = k(:)

    !check rates > 1
    if(maxval(k)>1d0) then
       print *,"ERROR: found rate >1d0 in coolingZ_rates!"
       print *," Tgas =",Tgas
       do i=1,nZrate
          if(k(i)>1d0) print *,i,k(i)
       end do
       stop
    end if

    !check rates <0
    if(minval(k)<0d0) then
       print *,"ERROR: found rate <0d0 in coolingZ_rates!"
       print *," Tgas =",Tgas
       do i=1,nZrate
          if(k(i)<0d0) print *,i,k(i)
       end do
       stop
    end if

  end function coolingZ_rates

  !**********************
  function coolingZ_rate_tabs(inTgas)
    use krome_commons
    implicit none
    real*8::inTgas,Tgas,coolingZ_rate_tabs(nZrate),k(nZrate)
    integer::idx,j
    Tgas = inTgas

    idx = (log10(Tgas)-coolTab_logTlow) * inv_coolTab_idx + 1

    idx = max(idx,1)
    idx = min(idx,coolTab_n-1)

    do j=1,nZrate
       k(j) = (Tgas-coolTab_T(idx)) * inv_coolTab_T(idx) * &
            (coolTab(j,idx+1)-coolTab(j,idx)) + coolTab(j,idx)
       k(j) = max(k(j), 0d0)
    end do

    coolingZ_rate_tabs(:) = k(:)

  end function coolingZ_rate_tabs

  !**********************
  subroutine coolingZ_init_tabs()
    use krome_commons
    implicit none
    integer::j,jmax,idx
    real*8::Tgas,Tgasold

    jmax = coolTab_n !size of the cooling tables (number of saples)

    !note: change upper and lower limit for rate tables here
    coolTab_logTlow = log10(2d0)
    coolTab_logTup = log10(1d8)

    !pre compute this value since used jmax times
    inv_coolTab_idx = (jmax-1) / (coolTab_logTup-coolTab_logTlow)

    !loop over the jmax interpolation points
    do j=1,jmax
       !compute Tgas for the given point
       Tgas = 1d1**((j-1)*(coolTab_logTup-coolTab_logTlow) &
            /(jmax-1) + coolTab_logTlow)
       !produce cooling rates for the given Tgas
       coolTab(:,j) = coolingZ_rates(Tgas)
       !store Tgas into the array
       coolTab_T(j) = Tgas
       !save 1/dT since it is known
       if(j>1) inv_coolTab_T(j-1) = 1d0 / (Tgas-Tgasold)
       Tgasold = Tgas
    end do

  end subroutine coolingZ_init_tabs

  !*******************************
  !this subroutine solves a non linear system
  ! with the equations stored in fcn function
  ! and a dummy jacobian jcn
  subroutine nleq_wrap(x)
    use krome_user_commons
    integer,parameter::nmax=100 !problem size
    integer,parameter::liwk=nmax+50 !size integer workspace
    integer,parameter::lrwk=(nmax+13)*nmax+60 !real workspace
    integer,parameter::luprt=6 !logical unit verbose output
    integer::neq,iopt(50),ierr,niw,nrw,iwk(liwk),ptype,i
    real*8::x(:),xscal(nmax),rtol,rwk(lrwk),idamp,mdamp,xi(size(x)),minx
    real*8::store_invdvdz
    neq = size(x)
    niw = neq+50
    nrw = (neq+13)*neq+60

    ptype = 2 !initial problem type, 2=mildly non-linear
    rtol = 1d-5 !realtive tolerance
    xi(:) = x(:) !store initial guess
    idamp = 1d-4 !initial damp (when ptype>=4, else default)
    mdamp = 1d-8 !minimum damp (when ptype>=4, else default)
    ierr = 0

    !iterate until ierr==0 and non-negative solutions
    do
       if(ptype>50) then
          print *,"ERROR in nleq1: can't find a solution after attempt",ptype
          stop
       end if
       
       x(:) = xi(:) !restore initial guess

       !if damping error or negative solutions
       ! prepares initial guess with the thin case
       if(ptype>7.and.(ierr==3.or.ierr==0)) then
          rtol = 1d-5
          iwk(:) = 0
          iopt(:) = 0
          rwk(:) = 0d0
          xscal(:) = 0d0
          store_invdvdz = krome_invdvdz !store global variable
          krome_invdvdz = 0d0 !this sets beta to 1
#IFKROME_use_NLEQ
          call nleq1(neq,fcn,jcn,x(:),xscal(:),rtol,iopt,ierr,&
               liwk,iwk(:),lrwk,rwk(:))
#ENDIFKROME_use_NLEQ
          if(ierr.ne.0) then
             print *,"ERROR in nleq for thin approx",ierr
             stop
          end if
          krome_invdvdz = store_invdvdz !restore global variable
       end if
       xscal(:) = 0d0 !scaling factor
       rtol = 1d-5 !relative tolerance
       iwk(:) = 0 !default iwk
       iwk(31) = int(1e8) !max iterations
       iopt(:) = 0 !default iopt
       iopt(31) = min(ptype,4) !problem type
       rwk(:) = 0d0 !default rwk
       !reduce damps if damping error
       if(ptype>4.and.ierr==3) then
          idamp = idamp * 1d-1 !reduce idamp
          mdamp = mdamp * 1d-1 !reduce mdamp
       end if
       !if problem is extremely nonlinear use custom damps
       if(ptype>4) then
          rwk(21) = idamp !copy idamp to solver
          rwk(22) = mdamp !copy mdamp to solver
       end if

#IFKROME_use_NLEQ
       call nleq1(neq,fcn,jcn,x(:),xscal(:),rtol,iopt,ierr,&
            liwk,iwk(:),lrwk,rwk(:))
#ENDIFKROME_use_NLEQ

       !check for errors
       if(ierr.ne.0) then
          !print *,"error",ierr
          !problem with damping factor and/or problem type
          if(ierr==3) then
             ptype = ptype + 1 !change the problem type (non-linearity)
          elseif(ierr==5) then
             xi(:) = x(:)
          else
             !other type of error hence stop
             print *,"ERROR in nleq1, ierr:",ierr
             print *,"solutions found so far:"
             do i=1,size(x)
                print *,i,x(i)
             end do
             stop
          end if
       else
          !if succesful search for negative results
          minx = minval(x) !minimum value
          !if minimum value is positive OK
          if(minx.ge.0d0) then
             exit
          else
             !if negative values are small set to zero
             if(abs(minx)/maxval(x)<rtol) then
                do i=1,neq
                   x(i) = max(x(i),0d0)
                end do
                exit
             else
                !if large negative values increase non-linearity
                ptype = ptype + 1
             end if
          end if
       end if
    end do
  end subroutine nleq_wrap

  !***************************
  subroutine fcn(n,x,f,ierr)
    implicit none
    integer::n,ierr
    real*8::x(n),f(n)

#KROME_fcn_cases

  end subroutine fcn

  !**********************************
  !dummy jacobian for non linear equation solver
  subroutine jcn()

  end subroutine jcn

#KROME_coolingZ_functions

#ENDIFKROME

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

#IFKROME_useLAPACK
  !*********************************
  subroutine mydgesv(n,Ain,Bin, parent_name)
    !driver for LAPACK dgesv
    integer::n,info,i,ipiv(n)
    real*8,allocatable::tmp(:)
    real*8::A(n,n),B(n),Ain(:,:),Bin(:),suml,sumr,tmpn(n)
    character(len=*)::parent_name
    A(:,:) = Ain(1:n,1:n)
    B(:) = Bin(1:n)
    call dgesv(n,1,A,n,ipiv,B,n,info)
    Bin(1:n) = B(:)

    !write some info about the error and stop
    if(info > 0) then
       allocate(tmp(size(Bin)))
       print *,"ERROR: matrix exactly singular, U(i,i) where i=",info
       print *,' (called by "'//trim(parent_name)//'" function)'
       
       !dump the input matrix to a file
       open(97,file="ERROR_dump_dgesv.dat",status="replace")
       !dump size of the problem
       write(97,*) "size of the problem:",n
       write(97,*)

       !dump matrix A
       write(97,*) "Input matrix A line by line:"
       do i=1,size(Ain,1)
          tmp(:) = Ain(i,:)
          write(97,'(I5,999E17.8e3)') i,tmp(:)
       end do

       !dump matrix A
       write(97,*)
       write(97,*) "Workin matrix A line by line:"
       do i=1,n
          tmpn(:) = Ain(i,1:n)
          write(97,'(I5,999E17.8e3)') i,tmpn(:)
       end do

       !dump matrix B
       write(97,*)
       write(97,*) "Input/output vector B element by element"
       do i=1,n
          write(97,*) i, Bin(i),B(i)
       end do

       !dump info on matrix A rows
       write(97,*)
       write(97,*) "Info on matrix A rows"
       write(97,'(a5,99a17)') "idx","minval","maxval"
       do i=1,size(Ain,1)
           write(97,'(I5,999E17.8e3)') i, minval(Ain(i,:)), &
                maxval(Ain(i,:))
       end do

       !dump info on matrix sum left and right
       write(97,*)
       write(97,*) "Info on matrix A, sum left/right"
       write(97,'(a5,99a17)') "idx","left","right"
       suml = 0d0
       sumr = 0d0
       do i=1,size(Ain,1)
          if(i>1) suml = sum(Ain(i,:i-1))
          if(i<n) sumr = sum(Ain(i,i+1:))
          write(97,'(I5,999E17.8e3)') i, suml, sumr
       end do
       close(97)

       print *,"Input A and B dumped in ERROR_dump_dgesv.dat"

       stop
    end if
    
    !if error print some info and stop
    if(info<0) then
       print *,"ERROR: input error position ",info
       print *,' (called by "'//trim(parent_name)//'" function)'
       stop
    end if

  end subroutine mydgesv
#ENDIFKROME

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
#IFKROME_useCoolingH2
       cool_H2 = cooling_H2(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingH2GP
       cool_H2GP = cooling_H2GP(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingAtomic
       cool_atomic = cooling_atomic(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingHD
       cool_HD = cooling_HD(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingZ
       cool_Z = cooling_Z(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingdH
       cool_dH = cooling_dH(n(:),Tgas)
#ENDIFKROME
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

#IFKROME_useCoolingH2
    cool(1) = cooling_H2(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingH2GP
    cool(2) = cooling_H2GP(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingAtomic
    cool(3) = cooling_atomic(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingHD
    cool(4) = cooling_HD(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingZ
    cool(5) = cooling_Z(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingdH
    cool(6) = cooling_dH(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDust
    cool(7) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCompton
    cool(8) = cooling_compton(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCIE
    cool(9) = cooling_CIE(n(:), Tgas)
#ENDIFKROME

    write(nfile,'(99E14.5e3)') Tgas,sum(cool),cool(:) 

  end subroutine dump_cool

end module KROME_cooling
