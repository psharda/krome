module KROME_cooling
#KROME_header
  integer,parameter::coolTab_n=int(1e3)
#KROME_nZrate
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx

contains

  !*******************************
  function cooling(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::cooling,cools(10)
    
    !returns cooling in erg/cm3/s
    cools(:) = 0.d0

#IFKROME_useCoolingH2
    cools(1) = cooling_H2(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingH2GP
    cools(2) = cooling_H2GP(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingAtomic
    cools(3) = cooling_Atomic(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingHD
    cools(4) = cooling_HD(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingZ
    cools(5) = cooling_Z(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingdH
    cools(6) = cooling_dH(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingDust
    cools(7) = cooling_dust(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCompton
    cools(8) = cooling_compton(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingCIE
    cools(9) = cooling_CIE(n(:), Tgas)
#ENDIFKROME

#IFKROME_useCoolingContinuum
    cools(10) = cooling_continuum(n(:), Tgas)
#ENDIFKROME

    cooling = sum(cools)

    !remove the comment below to write cooling contributions to fort.44
    !write(44,'(99E17.8e3)') sum(n(1:nmols)),Tgas,cools(:)

    !gnuplot command (n=100, and m=1 for density or m=2 for temperature) 
    !plot 'fort.44' u m:3 every n w l t "H2",\
    ! '' u m:4 every n w l t "H2GP",\
    ! '' u m:5 every n w l t "atomic",\
    ! '' u m:6 every n w l t "HD",\
    ! '' u m:7 every n w l t "Z",\
    ! '' u m:8 every n w l t "enthalpy",\
    ! '' u m:9 every n w l t "dust",\
    ! '' u m:10 every n w l t "compton",\
    ! '' u m:11 every n w l t "CIE",\
    ! '' u m:12 every n w l t "Cont"

  end function cooling


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

  !*****************************
  function coolingChem(n,Tgas)
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
  function cooling_CIE(n, Tgas)
    !CIE cooling: fit from Ripamponti&Abel2004 (RA04) data
    ! The fit is valid from 100K-1e6K.
    ! Original data are from 400K to 7000K.
    ! We extrapolated data under 400K and fitted from 100K to 10**2.95 K.
    ! Data from 10**2.95 K to 1e5K are fitted analogously.
    ! Above 1e5 we employ a cubic extrapolation.
    use krome_commons
    use krome_constants
    real*8::cooling_CIE,n(:),Tgas
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
    Tgas = max(Tgas, 2.73d0) 

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


#IFKROME_useCoolingCompton
  !*******************************
  function cooling_compton(n, Tgas)
    !compton cooling erg/cm3/s from Cen1992
    use krome_user_commons
    use krome_commons
    real*8::cooling_compton,n(:),Tgas
    
    !note that redhsift is defined in krome_user_commons and 
    ! must be provided by the user
    cooling_compton = 5.65d-36 * (1.d0 + redshift)**4 &
         * (Tgas - 2.73d0 * (1.d0 + redshift)) * n(idx_e) !erg/s/cm3
    
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
    real*8::logt32,logt33,logt34,logt35
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

    !//H2-H
    if(temp>1d1 .and. temp<=1d2) then
       cool = cool +1.d1**(-16.818342D0 +3.7383713D1*logt3 &
            +5.8145166D1*logt32 +4.8656103D1*logt33 &
            +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
       cool = cool +1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
            -1.1332860D1*logt32 -2.7850082D1*logt33 &
            -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
    elseif(temp>1.d3 .and. temp<=6.d3) then
       cool = cool +1d1**(-2.4311209D1 +4.6450521D0*logt3 &
            -3.7209846D0*logt32 +5.9369081D0*logt33 &
            -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
    end if

    !//H2-Hp
    if(temp>1.d1 .and. temp<=1.d4)  then
       cool = cool + 1d1**(-2.1716699D1 +1.3865783D0*logt3 &
            -0.37915285D0*logt32 +0.11453688D0*logt33 &
            -0.23214154D0*logt34 +0.058538864D0*logt35)*n(idx_Hj)
    end if

    !//H2-H2
    if(temp>1.d2 .and. temp<=6.d3) then
       cool = cool + 1d1**(-2.3962112D1 +2.09433740D0*logt3 &
            -.77151436D0*logt32 +.43693353D0*logt33 &
            -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2)
    end if

    !//H2-e
    if(temp>1d1 .and. temp<=2d2) then
       cool =  cool +1d1**(-3.4286155D1 -4.8537163D1*logt3 &
            -7.7121176D1*logt32 -5.1352459D1*logt33 &
            -1.5169150D1*logt34 -.98120322D0*logt35)*n(idx_e)
    elseif(temp>2d2 .and. temp<1d4)  then
       cool = cool + 1d1**(-2.2190316D1 +1.5728955D0*logt3 &
            -.213351D0*logt32 +.96149759D0*logt33 &
            -.91023195D0*logt34 +.13749749D0*logt35)*n(idx_e)
    end if

    !//H2-He
    if(temp>1.d1 .and. temp<=6.d3) then
       cool =  cool + 1d1**(-2.3689237d1 +2.1892372d0*logt3&
            -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
            +.19191375d0*logt35)*n(idx_He)
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
    real*8::temp,gaunt_factor,T5,cool,bms_ions


    gaunt_factor = 1.5d0
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

    !BREMSSTRAHLUNG: all ions
#KROME_brem_ions
    cool = cool+ 1.42d-27*gaunt_factor*sqrt(temp)&
         *bms_ions*n(idx_e)

    cooling_atomic = max(cool, 0.d0)  !erg/cm3/s
  end function cooling_Atomic
#ENDIFKROME

#IFKROME_useCoolingHD
  !HD COOLING LIPOVKA ET AL. MNRAS, 361, 850, (2005)
  !UNITS=erg/cm3/s
  !*******************************
  function cooling_HD(n, Tgas)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer, parameter::ns=4
    real*8::cooling_HD
    real*8::n(:),Tgas,logTgas,lognH
    real*8::c(0:ns,0:ns),logW,W,dd,lhj

    !default HD cooling value
    cooling_HD = 0.0d0 !erg/cm3/s

    !exit on low temperature
    if(Tgas<1d2) return

    !calculate density
    dd = n(idx_H) !sum(n(1:nmols))
    !exit if density is out of Lipovka bounds
    if(dd<1d0 .or. dd>1d8) return 

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
    do i=1,n-1
       if(x.le.xval(i)) then
          !linear fit
          flin = (yval(i+1) - yval(i)) / (xval(i+1) - xval(i)) * &
               (x - xval(i)) + yval(i)
          found = .true. !found flag
          exit
       end if
    end do
    if(.not.found) flin = yval(n)

  end function flin

  !***********************
  !metal cooling as in Maio et al. 2007
  ! loaded from data file 
  function cooling_Z(n,inTgas)
    !//coefficients written as gijCOOLANT_COLLIDER
    !//e.g. g10C_H2o means Carbon colliding with H2-ortho cooling through 
    !//transition 1->0
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
    implicit none
    real*8::inTgas, coolingZ_rates(nZrate),k(nZrate)
    real*8::Tgas
#KROME_coolingZ_declare_custom_vars

    Tgas = inTgas

#KROME_coolingZ_custom_vars

#KROME_coolingZ_rates

    coolingZ_rates(:) = k(:)

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
    coolTab_logTlow = log10(3d0)
    coolTab_logTup = log10(1d4)

    !pre compute this value since used jmax times
    inv_coolTab_idx = (jmax-1) / (coolTab_logTup-coolTab_logTlow)
    
    !loop over the jmax interpolation points
    do j=1,jmax
       !compute Tgas for the given point
       Tgas = 1d1**((j-1)*(coolTab_logTup-coolTab_logTlow) &
            /(jmax-1) + coolTab_logTlow)
       !produce cooling function for the given Tgas
       coolTab(:,j) = coolingZ_rates(Tgas)
       !store Tgas into the array
       coolTab_T(j) = Tgas
       !save 1/dT since it is known
       if(j>1) inv_coolTab_T(j-1) = 1d0 / (Tgas-Tgasold)
       Tgasold = Tgas
    end do

  end subroutine coolingZ_init_tabs

#KROME_coolingZ_functions


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

  !*********************************
  subroutine mydgesv(A,B)
    !driver for LAPACK dgesv
    real*8::A(:,:),B(:)
    real*8,allocatable::tmp(:)
    integer::n,info,i
    integer,allocatable::ipiv(:)
    n=size(B)
    allocate(tmp(n))
    allocate(ipiv(n))
    call dgesv(n,1,A,n,ipiv,B,n,info)
    if(info .ne. 0) then
       print *,"ERROR: problem with dgesv, info: ",info
       do i=1,n
          tmp(:)=A(i,:)
          write(*,*) tmp(:)
       end do
       write(*,*) "B:",B(:)
       stop
    end if
    deallocate(tmp)
    deallocate(ipiv)
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
