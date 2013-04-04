module KROME_cooling
contains

#KROME_header


  !*******************************
  function cooling(n, Tgas)
    implicit none
    real*8::n(:), Tgas
    real*8::cooling 
    !total cooling erg/cm3/s
    cooling = 0.d0

#IFKROME_useCoolingH2
    cooling = cooling + cooling_H2(n(:), Tgas)
#ENDIFKROME
#IFKROME_useCoolingH2GP
    cooling = cooling + cooling_H2GP(n(:), Tgas)
#ENDIFKROME
#IFKROME_useCoolingCEN
    cooling = cooling + cooling_CEN(n(:), Tgas)
#ENDIFKROME
#IFKROME_useCoolingHD
    cooling = cooling + cooling_HD(n(:), Tgas)
#ENDIFKROME
#IFKROME_useCoolingZ
    cooling = cooling + cooling_Z(n(:), Tgas)
#ENDIFKROME

  end function cooling

  !*******************************
  function heating(n, Tgas)
    implicit none
    real*8::n(:), Tgas
    real*8::heating 
    !total heating erg/cm3/s
    heating = 0.d0
#IFKROME_useHeatingA
    heating = heating + heatingA(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingCompress
    heating = heating + heat_compress(n(:), Tgas)
#ENDIFKROME

#IFKROME_useHeatingPhoto
    heating = heating + photo_heating(n(:))
#ENDIFKROME
    
  end function heating

#IFKROME_useCoolingH2GP
  !H2 COOLING GALLI&PALLA 1998 
  !VALIDITY RANGE 10<=T(K)<=1d4
  !*******************************
  function cooling_H2GP(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::n(:),Tgas, tm, logT
    real*8::cooling_H2GP,T3
    real*8::LDL,HDLR,HDLV,HDL,fact

    tm = max(Tgas, 13.0d0)    ! no cooling below 13 Kelvin ...
    tm = min(Tgas, 1.d5)      ! fixes numerics
    logT = log10(tm)
    T3 = tm / 1.d3

    !low density limit in erg/s
    LDL = 1.d1**(-103.d0+97.59d0*logT-48.05d0*logT**2&
         +10.8d0*logT**3-0.9032d0*logT**4)*n(idx_H)

    !high density limit
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
         3.e-24*exp(-0.51/t3)) !erg/s
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
    HDL  = HDLR + HDLV !erg/s

    !TO AVOID DIVISION BY ZERO
    if(LDL.eq.0.d0)then
       fact = 0.0d0
    else
       fact = HDL/LDL !dimensionless
    endif

    cooling_H2GP = HDL*n(idx_H2)/(1.d0+(fact))!erg/cm3/s

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
    integer::i
    character*16::names(nspec)
    temp = max(Tgas, 1d1)

    T3 = temp / 1.d3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0.d0

    !//H2-H
    if(temp>1d1 .and. temp<=1d2) then
       cool = cool +1.d1**(-16.818342D0 +3.7383713D1*logt3 &
            +5.8145166D1*logt3**2 +4.8656103D1*logt3**3 &
            +2.0159831D1*logt3**4 +3.8479610D0*logt3**5 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
       cool = cool +1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
            -1.1332860D1*logt3**2 -2.7850082D1*logt3**3 &
            -2.1328264D1*logt3**4 -4.2519023D0*logt3**5)*n(idx_H)
    elseif(temp>1.d3 .and. temp<=6.d3) then
       cool = cool +1d1**(-2.4311209D1 +4.6450521D0*logt3 &
            -3.7209846D0*logt3**2 +5.9369081D0*logt3**3 &
            -5.5108049D0*logt3**4 +1.5538288D0*logt3**5)*n(idx_H)
    end if

    !//H2-Hp
    if(temp>1.d1 .and. temp<=1.d4)  then
       cool = cool + 1d1**(-2.1716699D1 +1.3865783D0*logt3 &
            -0.37915285D0*logt3**2 +0.11453688D0*logt3**3 &
            -0.23214154D0*logt3**4 +0.058538864D0*logt3**5)*n(idx_Hj)
    end if

    !//H2-H2
    if(temp>1.d2 .or. temp<=6.d3) then
       cool = cool + 1d1**(-2.3962112D1 +2.09433740D0*logt3 &
            -.77151436D0*logt3**2 +.43693353D0*logt3**3 &
            -.14913216D0*logt3**4 -.033638326D0*logt3**5)*n(idx_H2)
    end if

    !//H2-e
    if(temp>1d1 .and. temp<=2d2) then
       cool =  cool +1d1**(-3.4286155D1 -4.8537163D1*logt3 &
            -7.7121176D1*logt3**2 -5.1352459D1*logt3**3 &
            -1.5169150D1*logt3**4 -.98120322D0*logt3**5)*n(idx_e)
    elseif(temp>2d2 .and. temp<1d4)  then
       cool = cool + 1d1**(-2.2190316D1 +1.5728955D0*logt3 &
            -.213351D0*logt3**2 +.96149759D0*logt3**3 &
            -.91023195D0*logt3**4 +.13749749D0*logt3**5)*n(idx_e)
    end if

    !//H2-He
    if(temp>1.d1 .and. temp<=6.d3) then
       cool =  cool + 1d1**(-2.3689237d1 +2.1892372d0*logt3&
            -.81520438d0*logt3**2 +.29036281d0*logt3**3 -.16596184d0*logt3**4 &
            +.19191375d0*logt3**5)*n(idx_He)
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

    !high density limit from HM79, GP98 
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL 
    HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
         3.e-24*exp(-0.51/t3))
    HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3))
    HDL  = HDLR + HDLV

    LDL = cool !erg/s

    !TO AVOID DIVISION BY ZERO
    if(LDL.eq.0.d0)then
       fact = 0.0d0
    else
       fact = HDL/LDL 
    endif

    cooling_H2 = HDL*n(idx_H2)/(1.d0+(fact)) !erg/cm3/s

  end function cooling_H2
#ENDIFKROME


#IFKROME_useCoolingCEN
  !CEN COOLING ApJS, 78, 341, 1992
  !UNITS = erg/s/cm3
  !*******************************
  function cooling_CEN(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::Tgas,cooling_CEN,n(:)
    real*8::temp,gaunt_factor,T5,cool


    gaunt_factor = 1.5d0
    temp = max(Tgas,10.d0) !K
    T5 = temp/1.d5 !K
    cool = 0.d0

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
    cool = cool+ 1.42d-27*gaunt_factor*sqrt(temp)&
         *(n(idx_Hj)+n(idx_H2j)+n(idx_Hej))*n(idx_e)

    cooling_CEN = max(cool, 0.d0)  !erg/cm3/s
  end function cooling_CEN
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
    cooling_HD = 0.0d0

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
    logW = 0.0d0
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
  !metal cooling (as in Maio et al. 2007)
  function cooling_Z(n,inTgas)
    !//coefficients written as gijCOOLANT_COLLIDER
    !//e.g. g10C_H2o means Carbon colliding with H2-ortho cooling through 
    !//transition 1->0
    use krome_commons
    use krome_constants
    implicit none
    real*8::cooling_Z,n(:)
    real*8::inTgas,Tgas,cool,lnT,kb
    real*8::nH2p,nH2o
    integer::i

    real*8::g10C_H2o,g20C_H2o,g21C_H2o
    real*8::g10C_H2p,g20C_H2p,g21C_H2p
    real*8::g10C_H,g20C_H,g21C_H
    real*8::g10C_Hp,g20C_Hp,g21C_Hp
    real*8::g10C_e,g20C_e,g21C_e
    real*8::g10Si_H, g20Si_H, g21Si_H
    real*8::g10Si_Hp, g20Si_Hp, g21Si_Hp
    real*8::g10Fe_H,g20Fe_H,g21Fe_H
    real*8::g10Fe_e,g20Fe_e,g21Fe_e
    real*8::g30Fe_e,g40Fe_e,g43Fe_e
    real*8::g10O_H2o,g20O_H2o,g21O_H2o
    real*8::g10O_H2p,g20O_H2p,g21O_H2p
    real*8::g10O_H,g20O_H,g21O_H
    real*8::g10O_Hp,g20O_Hp,g21O_Hp
    real*8::g10O_e,g20O_e,g21O_e
    real*8::g10Cp_H,g10Cp_e
    real*8::g10Op_e,g20Op_e,g21Op_e
    real*8::g10Sip_e,g10Sip_H
    real*8::g10Fep_H,g21Fep_H,g32Fep_H
    real*8::g43Fep_H,g20Fep_H,g30Fep_H
    real*8::g40Fep_H,g31Fep_H,g41Fep_H,g42Fep_H
    real*8::g10Fep_e,g21Fep_e,g32Fep_e
    real*8::g43Fep_e,g20Fep_e,g30Fep_e
    real*8::g40Fep_e,g31Fep_e,g41Fep_e,g42Fep_e

    real*8::g01C_H2o,g02C_H2o,g12C_H2o
    real*8::g01C_H2p,g02C_H2p,g12C_H2p
    real*8::g01C_H  ,g02C_H   ,g12C_H
    real*8::g01C_Hp ,g02C_Hp  ,g12C_Hp
    real*8::g01C_e  ,g02C_e   ,g12C_e
    real*8::g01Si_H, g02Si_H, g12Si_H
    real*8::g01Si_Hp, g02Si_Hp, g12Si_Hp
    real*8::g01Fe_H,g02Fe_H,g12Fe_H
    real*8::g01Fe_e,g02Fe_e,g12Fe_e
    real*8::g03Fe_e,g04Fe_e,g34Fe_e
    real*8::g01O_H2o,g02O_H2o,g12O_H2o
    real*8::g01O_H2p,g02O_H2p,g12O_H2p
    real*8::g01O_H  ,g02O_H  ,g12O_H
    real*8::g01O_Hp ,g02O_Hp ,g12O_Hp
    real*8::g01O_e  ,g02O_e  ,g12O_e
    real*8::g01Cp_H ,g01Cp_e
    real*8::g01Op_e ,g02Op_e ,g12Op_e
    real*8::g01Sip_e,g01Sip_H
    real*8::g01Fep_H,g12Fep_H,g23Fep_H
    real*8::g34Fep_H,g02Fep_H,g03Fep_H
    real*8::g04Fep_H,g13Fep_H,g14Fep_H,g24Fep_H
    real*8::g01Fep_e,g12Fep_e,g23Fep_e
    real*8::g34Fep_e,g02Fep_e,g03Fep_e
    real*8::g04Fep_e,g13Fep_e,g14Fep_e,g24Fep_e

    real*8::C_g10to01,C_g20to02,C_g21to12
    real*8::Si_g10to01,Si_g20to02,Si_g21to12
    real*8::Fe_g10to01,Fe_g20to02,Fe_g21to12
    real*8::O_g10to01,O_g21to12,O_g20to02
    real*8::Cp_g10to01
    real*8::Sip_g10to01

    real*8::M01C,M12C,M02C
    real*8::M10C,M21C,M20C
    real*8::M01Si,M12Si,M02Si
    real*8::M10Si,M21Si,M20Si
    real*8::M01Fe,M12Fe,M02Fe,M03Fe,M04Fe,M34Fe
    real*8::M10Fe,M21Fe,M20Fe,M30Fe,M40Fe,M43Fe
    real*8::M01O,M12O,M02O
    real*8::M10O,M21O,M20O
    real*8::M01Cp,M10Cp
    real*8::M01Op,M02Op,M12Op
    real*8::M10Op,M20Op,M21Op
    real*8::M01Sip,M10Sip
    real*8::M01Fep,M12Fep,M23Fep
    real*8::M34Fep,M02Fep,M03Fep
    real*8::M04Fep,M13Fep,M14Fep,M24Fep
    real*8::M10Fep,M21Fep,M32Fep
    real*8::M43Fep,M20Fep,M30Fep
    real*8::M40Fep,M31Fep,M41Fep,M42Fep

    real*8::AC(3,3),BC(3)
    real*8::ASi(3,3),BSi(3)
    real*8::AFe(5,5),BFe(5)
    real*8::AO(3,3),BO(3)
    real*8::ACp(2,2),BCp(2)
    real*8::AOp(3,3),BOp(3)
    real*8::ASip(2,2),BSip(2)
    real*8::AFep(5,5),BFep(5)

    real*8::tot_metals,invTgas,ratio_para_ortho

    if(minval(n)<0.d0)then
       do i=1,size(n)
          if(n(i)<0.d0) n(i)=0.d0
       end do
    end if

    tot_metals = n(idx_C) + n(idx_Cj) + n(idx_Si) + n(idx_Sij) + n(idx_Fe)&
         + n(idx_Fej) + n(idx_O) + n(idx_Oj)
    ratio_para_ortho = 1.d0/3.d0
    kb = boltzmann_erg
    Tgas = max(inTgas,1.d1)
    invTgas = 1.d0/Tgas
    lnT = log(Tgas)
    cool = 0.d0
    if(Tgas.le.1d4 .and. tot_metals.ge.1.d-20) then
       !(1)------------------------------------------------------------------
       !//C-H2_ortho
       g10C_H2o = 8.7d-11 -6.6d-11*exp(-Tgas/218.3)+6.6d-11*exp(-2.*Tgas/218.3)
       g20C_H2o = 1.2d-10 -6.1d-11*exp(-Tgas/387.3)
       g21C_H2o = 2.9d-10 -1.9d-10*exp(-Tgas/348.9)

       !//C-H2_para
       g10C_H2p = 7.9D-11 -8.7D-11*EXP(-Tgas/126.4) +1.3D-10*EXP(-2.*Tgas/126.4)
       g20C_H2p = 1.1D-10 -8.6D-11*EXP(-Tgas/233.) +8.7D-11*EXP(-2.*Tgas/233.)
       g21C_H2P = 2.7D-10 -2.6D-10*EXP(-Tgas/250.7) +1.8D-10*EXP(-2.*Tgas/250.7)

       !//C-H
       g10C_H = 1.6D-10*(Tgas/100.)**(.14)
       g20C_H = 9.2D-11*(Tgas/100.)**(.26)
       g21C_H = 2.9D-10*(Tgas/100.)**(.26)

       !//C-Hp
       g10C_Hp = (9.6D-11 -1.8D-14*Tgas +1.9D-18*Tgas**2) *Tgas**(.45)
       IF(Tgas > 5d3)  g10C_Hp = 8.9D-10*Tgas**(.117)
       g20C_Hp = (3.1D-12 -6.D-16*Tgas +3.9d-20*Tgas**2) *Tgas
       IF(Tgas > 5d3)  g20C_Hp = 2.3D-9*Tgas**(.0965)
       g21C_Hp = (1.D-10 -2.2D-14*Tgas +1.7D-18*Tgas**2) *Tgas**(.70)
       IF(Tgas > 5.d3)  g21C_Hp = 9.2D-9*Tgas**(.0535)

       !//C-e
       g10C_e = 2.88D-6*Tgas**(-.5)*EXP(-9.25141 -7.73782D-1*lnT &
            +3.61184D-1*lnT**2 &
            -1.50892D-2*lnT**3 -6.56325D-4*lnT**4)
       IF(Tgas > 1D3)  g10C_e = 2.88D-6*Tgas**(-.5) *EXP(4.446D2 &
            -2.27913D2*lnT +4.2595D1*lnT**2 -3.4762*lnT**3 +1.0508D-1*lnT**4)
       g20C_e = 1.73D-6*Tgas**(-.5)*EXP(-7.69735 -1.30743*lnT +.697638*lnT**2 &
            -.111338*lnT**3 +.705277D-2*lnT**4)
       IF(Tgas > 1D3)  g20C_e = 1.73D-6*Tgas**(-.5)*EXP(3.50609D2 &
            -1.87474D2*lnT +3.61803D1*lnT**2 -3.03283*lnT**3 +9.38138D-2*lnT**4)
       g21C_e = 1.73D-6*Tgas**(-.5)*EXP(-7.4387 -.57443*lnT +.358264*lnT**2 &
            -4.18166D-2*lnT**3 +2.35272D-3*lnT**4)
       IF(Tgas > 1D3)  g21C_e = 1.73D-6*Tgas**(-.5)*EXP(3.86186D2 &
            -2.02192D2*lnT +3.85049D1*lnT**2 -3.19268*lnT**3 +9.78573D-2*lnT**4)

       !//Si-H
       g10Si_H = 3.5D-10*(Tgas/1D2)**(-.03)
       g20Si_H = 1.7D-11*(Tgas/1D2)**(.17)
       g21Si_H = 5.D-10*(Tgas/1.D2)**(.17)

       !//Si-H+
       g10Si_Hp = 7.2D-9
       g20Si_Hp = 7.2D-9
       g21Si_Hp = 2.2D-8

       !//Fe-H
       g10Fe_H = 8.D-10*(Tgas/1.D2)**(.17)
       g20Fe_H = 6.9D-10*(Tgas/1.D2)**(.17)
       g21Fe_H = 5.3D-10*(Tgas/1.D2)**(.17)

       !//Fe-e
       g10Fe_e = 1.2D-7
       g20Fe_e = 1.2D-7
       g21Fe_e = 9.3D-8
       g30Fe_e = 2.D-7*(Tgas/1.D4)**(.57)
       IF(Tgas > 1.D4) g30Fe_e = 2.D-7*(Tgas/1.D4)**(-.13)
       g40Fe_e = 1.D-7*(Tgas/1.D4)**(.57)
       IF(Tgas > 1.D4) g40Fe_e = 1.D-7
       g43Fe_e = 1.5D-7

       !//O-H2o
       g10O_H2o = 2.7D-11*Tgas**(.362)
       g20O_H2o = 5.49D-11*Tgas**(.317)
       g21O_H2o = 2.74D-14*Tgas**(1.06)

       !//O-H2p
       g10O_H2p = 3.46D-11*Tgas**(.316)
       g20O_H2p = 7.07D-11*Tgas**(.268)
       g21O_H2p = 3.33D-15*Tgas**(1.36)

       !//O-H
       g10O_H = 9.2D-11*(Tgas/100.)**(.67)
       g20O_H = 4.3D-11*(Tgas/100.)**(.80)
       g21O_H = 1.1D-10*(Tgas/100.)**(.44)

       !//O-Hp
       g10O_Hp = 6.38D-11*Tgas**(.4)
       IF(Tgas > 194.)  g10O_Hp = 7.75D-12*Tgas**(.8)
       IF(Tgas > 3686.)  g10O_Hp = 2.65D-10*Tgas**(.37)
       g20O_Hp = 6.1D-13*Tgas**(1.1)
       IF(Tgas > 511.) g20O_Hp = 2.12D-12*Tgas**(.9)
       IF(Tgas > 7510.) g20O_Hp = 4.49D-10*Tgas**(.3)
       g21O_Hp = 2.03D-11*Tgas**(.56)
       IF(Tgas > 2090.) g21O_Hp = 3.43D-10*Tgas**(.19)

       !//O-e
       g10O_e = 5.12D-10*Tgas**(-.075)
       g20O_e = 4.86D-10*Tgas**(-.026)
       g21O_e = 1.08D-14*Tgas**(.926)

       !//Cp-H
       g10Cp_H = 8D-10*(Tgas/100.)**(.07)
       !//Cp-e
       g10Cp_e = 2.8D-7*(Tgas/100.)**(-.5)

       !//Op-e
       g10Op_e = 1.3D-8*(Tgas/1.D4)**(-.5)
       g20Op_e = 1.3D-8*(Tgas/1.D4)**(-.5)
       g21Op_e = 2.5D-8*(Tgas/1.D4)**(-.5)

       !//Sip-e
       g10Sip_e = 1.2D-6*(Tgas/100.)**(-.5)
       !//Sip-H
       g10Sip_H = 4.95D-10*(Tgas/100.)**(.24)

       !//Fep-H
       g10Fep_H = 9.5D-10
       g21Fep_H = 4.7D-10
       g32Fep_H = 5.D-10
       g43Fep_H = 5.D-10
       g20Fep_H = 5.7D-10
       g30Fep_H = 5.D-10
       g40Fep_H = 5.D-10
       g31Fep_H = 5.D-10
       g41Fep_H = 5.D-10
       g42Fep_H = 5.D-10

       !//Fep-e
       g10Fep_e = 1.8D-6*(Tgas/100.)**(-.5)
       g21Fep_e = 8.7D-7*(Tgas/100.)**(-.5)
       g32Fep_e = 1.D-5*Tgas**(-.5)
       g43Fep_e = 1.D-5*Tgas**(-.5)
       g20Fep_e = 1.8D-6*(Tgas/100.)**(-.5)
       g30Fep_e = 1.D-5*Tgas**(-.5)
       g40Fep_e = 1.D-5*Tgas**(-.5)
       g31Fep_e = 1.D-5*Tgas**(-.5)
       g41Fep_e = 1.D-5*Tgas**(-.5)
       g42Fep_e = 1.D-5*Tgas**(-.5)


       !(2)-----------------------------------------------------------------
       !C:10->01
       C_g10to01 = 3./1.*EXP(-24.*invTgas)
       g01C_H2o = g10C_H2o*C_g10to01
       g01C_H2p = g10C_H2p*C_g10to01
       g01C_H = g10C_H*C_g10to01
       g01C_Hp = g10C_Hp*C_g10to01
       g01C_e = g10C_e*C_g10to01

       !//C:20->02
       C_g20to02 = 5./1.*EXP(-63.*invTgas)
       g02C_H2o = g20C_H2o*C_g20to02
       g02C_H2p = g20C_H2p*C_g20to02
       g02C_H = g20C_H*C_g20to02
       g02C_Hp = g20C_Hp*C_g20to02
       g02C_e = g20C_e*C_g20to02

       !//C:21->12
       C_g21to12 = 5./3.*EXP(-39.*invTgas)
       g12C_H2o = g21C_H2o*C_g21to12
       g12C_H2p = g21C_H2p*C_g21to12
       g12C_H = g21C_H*C_g21to12
       g12C_Hp = g21C_Hp*C_g21to12
       g12C_e = g21C_e*C_g21to12

       !//Si:10->01
       Si_g10to01 = 3./1.*EXP(-110.*invTgas)
       g01Si_H = g10Si_H*Si_g10to01
       g01Si_Hp = g10Si_Hp*Si_g10to01

       !//Si:20->02
       Si_g20to02 = 5./1.*EXP(-320.*invTgas)
       g02Si_H = g20Si_H*Si_g20to02
       g02Si_Hp = g20Si_Hp*Si_g20to02

       !//Si:21->12
       Si_g21to12 = 5./3.*EXP(-210.*invTgas)
       g12Si_H = g21Si_H*Si_g21to12
       g12Si_Hp = g21Si_Hp*Si_g21to12

       !//Fe:10->01
       Fe_g10to01 = 9./7.*EXP(-600.*invTgas)
       g01Fe_H = g10Fe_H*Fe_g10to01
       g01Fe_e = g10Fe_e*Fe_g10to01

       !//Fe:20->02
       Fe_g20to02 = 9./5.*EXP(-1000.*invTgas)
       g02Fe_H = g20Fe_H*Fe_g20to02
       g02Fe_e = g20Fe_e*Fe_g20to02

       !//Fe:21->12
       Fe_g21to12 = 7./5.*EXP(-420.*invTgas)
       g12Fe_H = g21Fe_H*Fe_g21to12
       g12Fe_e = g21Fe_H*Fe_g21to12

       !//Fe:ij->jk
       g03Fe_e = g30Fe_e*9./11.*EXP(-9900.*invTgas)
       g04Fe_e = g40Fe_e*9./8.*EXP(-11000.*invTgas)
       g34Fe_e = g43Fe_e*11./8.*EXP(-640.*invTgas)

       !//O:10->01
       O_g10to01 = 3./5.*EXP(-230.*invTgas)
       g01O_H2o = g10O_H2o*O_g10to01
       g01O_H2p = g10O_H2p*O_g10to01
       g01O_H = g10O_H*O_g10to01
       g01O_Hp = g10O_Hp*O_g10to01
       g01O_e = g10O_e*O_g10to01

       !//O:20->02
       O_g20to02 = 1./5.*EXP(-330.*invTgas)
       g02O_H2o = g20O_H2o*O_g20to02
       g02O_H2p = g20O_H2p*O_g20to02
       g02O_H = g20O_H*O_g20to02
       g02O_Hp = g20O_Hp*O_g20to02
       g02O_e = g20O_e*O_g20to02

       !//O:21->12
       O_g21to12 = 1./3.*EXP(-98.*invTgas)
       g12O_H2o = g21O_H2o*O_g21to12
       g12O_H2p = g21O_H2p*O_g21to12
       g12O_H = g21O_H*O_g21to12
       g12O_Hp = g21O_Hp*O_g21to12
       g12O_e = g21O_e*O_g21to12

       !//Cp:10->01
       Cp_g10to01 = 4./2.*EXP(-91.2*invTgas)
       g01Cp_e = g10Cp_e*Cp_g10to01
       g01Cp_H = g10Cp_H*Cp_g10to01

       !//Op:ij->ji
       g01Op_e = g10Op_e*11./7.*EXP(-39000.*invTgas)
       g02Op_e = g20Op_e*7./7.*EXP(-39000.*invTgas)
       g12Op_e = g21Op_e*7./11.*EXP(-30.*invTgas)

       !//Sip:10->01
       Sip_g10to01 = 2./4.*EXP(-410.*invTgas)
       g01Sip_e = g10Sip_e*Sip_g10to01
       g01Sip_H = g10Sip_H*Sip_g10to01

       !//Fep-H ij->ji
       g01Fep_H = g10Fep_H*10./8.*EXP(-553.58*invTgas)
       g12Fep_H = g21Fep_H*6./8.*EXP(-407.01*invTgas)
       g23Fep_H = g32Fep_H*6./4.*EXP(-280.57*invTgas)
       g34Fep_H = g43Fep_H*4./2.*EXP(-164.6*invTgas)
       g02Fep_H = g20Fep_H*10./6.*EXP(-960.59*invTgas)
       g03Fep_H = g30Fep_H*10./4.*EXP(-1241.16*invTgas)
       g04Fep_H = g40Fep_H*10./2.*EXP(-1405.76*invTgas)
       g13Fep_H = g31Fep_H*8./4.*EXP(-687.58*invTgas)
       g14Fep_H = g41Fep_H*8./2.*EXP(-852.18*invTgas)
       g24Fep_H = g42Fep_H*8./6.*EXP(-445.17*invTgas)

       !//Fep-e ij->ji
       g01Fep_e = g10Fep_e*10./8.*EXP(-553.58*invTgas)
       g12Fep_e = g21Fep_e*6./8.*EXP(-407.01*invTgas)
       g23Fep_e = g32Fep_e*6./4.*EXP(-280.57*invTgas)
       g34Fep_e = g43Fep_e*4./2.*EXP(-164.6*invTgas)
       g02Fep_e = g20Fep_e*10./6.*EXP(-960.59*invTgas)
       g03Fep_e = g30Fep_e*10./4.*EXP(-1241.16*invTgas)
       g04Fep_e = g40Fep_e*10./2.*EXP(-1405.76*invTgas)
       g13Fep_e = g31Fep_e*8./4.*EXP(-687.58*invTgas)
       g14Fep_e = g41Fep_e*8./2.*EXP(-852.18*invTgas)
       g24Fep_e = g42Fep_e*8./6.*EXP(-445.17*invTgas)


       !(3)-------------------------------------------------------------------
       !//calculates Mij=sum_k(n_k*g_ij**k)
       nH2p=n(idx_H2)/(ratio_para_ortho+1.) !para-H2
       nH2o=nH2p*ratio_para_ortho !orto-H2
       !//Carbon
       M01C = g01C_H2o*nH2o +g01C_H2p*nH2p +g01C_H*n(idx_H) +g01C_Hp*n(idx_Hj) &
            +g01C_e*n(idx_e)
       M12C = g12C_H2o*nH2o +g12C_H2p*nH2p +g12C_H*n(idx_H) +g12C_Hp*n(idx_Hj) &
            +g12C_e*n(idx_e)
       M02C = g02C_H2o*nH2o +g02C_H2p*nH2p +g02C_H*n(idx_H) +g02C_Hp*n(idx_Hj) &
            +g02C_e*n(idx_e)

       !//Si
       M01Si = g01Si_H*n(idx_H) +g01Si_Hp*n(idx_Hj)
       M12Si = g12Si_H*n(idx_H) +g12Si_Hp*n(idx_Hj)
       M02Si = g02Si_H*n(idx_H) +g02Si_Hp*n(idx_Hj)


       !//Fe
       M01Fe = g01Fe_H*n(idx_H) +g01Fe_e*n(idx_e)
       M12Fe = g12Fe_H*n(idx_H) +g12Fe_e*n(idx_e)
       M02Fe = g02Fe_H*n(idx_H) +g02Fe_e*n(idx_e)
       M03Fe = g03Fe_e*n(idx_e)
       M04Fe = g04Fe_e*n(idx_e)
       M34Fe = g34Fe_e*n(idx_e)

       !//O
       M01O = g01O_H2o*nH2o +g01O_H2p*nH2p +g01O_H*n(idx_H) +g01O_Hp*n(idx_Hj) &
            +g01O_e*n(idx_e)
       M12O = g12O_H2o*nH2o +g12O_H2p*nH2p +g12O_H*n(idx_H) +g12O_Hp*n(idx_Hj) &
            +g12O_e*n(idx_e)
       M02O = g02O_H2o*nH2o +g02O_H2p*nH2p +g02O_H*n(idx_H) +g02O_Hp*n(idx_Hj) &
            +g02O_e*n(idx_e)
       !//Cp
       M01Cp= g01Cp_H*n(idx_H) +g01Cp_e*n(idx_e)
       !//Op
       M01Op= g01Op_e*n(idx_e)
       M02Op= g02Op_e*n(idx_e)
       M12Op= g12Op_e*n(idx_e)
       !//Sip
       M01Sip= g01Sip_e*n(idx_e) +g01Sip_H*n(idx_H)
       !//Fep
       M01Fep= g01Fep_e*n(idx_e) +g01Fep_H*n(idx_H)
       M12Fep= g12Fep_e*n(idx_e) +g12Fep_H*n(idx_H)
       M23Fep= g23Fep_e*n(idx_e) +g23Fep_H*n(idx_H)
       M34Fep= g34Fep_e*n(idx_e) +g34Fep_H*n(idx_H)
       M02Fep= g02Fep_e*n(idx_e) +g02Fep_H*n(idx_H)
       M03Fep= g03Fep_e*n(idx_e) +g03Fep_H*n(idx_H)
       M04Fep= g04Fep_e*n(idx_e) +g04Fep_H*n(idx_H)
       M13Fep= g13Fep_e*n(idx_e) +g13Fep_H*n(idx_H)
       M14Fep= g14Fep_e*n(idx_e) +g14Fep_H*n(idx_H)
       M24Fep= g24Fep_e*n(idx_e) +g24Fep_H*n(idx_H)

       !//calculates Mji=sum_k(n_k*g_ji**k)+Ajk
       !//Carbon
       M10C = g10C_H2o*nH2o +g10C_H2p*nH2p +g10C_H*n(idx_H) +g10C_Hp*n(idx_Hj) &
            +g10C_e*n(idx_e) +7.9D-8
       M21C = g21C_H2o*nH2o +g21C_H2p*nH2p +g21C_H*n(idx_H) +g21C_Hp*n(idx_Hj) &
            +g21C_e*n(idx_e) +2.1D-14
       M20C = g20C_H2o*nH2o +g20C_H2p*nH2p +g20C_H*n(idx_H) +g20C_Hp*n(idx_Hj) &
            +g20C_e*n(idx_e) +2.7D-7

       !//Si
       M10Si = g10Si_H*n(idx_H) +g10Si_Hp*n(idx_Hj) +8.4D-6
       M21Si = g21Si_H*n(idx_H) +g21Si_Hp*n(idx_Hj) +2.4D-10
       M20Si = g20Si_H*n(idx_H) +g20Si_Hp*n(idx_Hj) +4.2D-5

       !//Fe
       M10Fe = g10Fe_H*n(idx_H) +g10Fe_e*n(idx_e) +2.5D-3
       M21Fe = g21Fe_H*n(idx_H) +g21Fe_e*n(idx_e) +1.6D-3
       M20Fe = g20Fe_H*n(idx_H) +g20Fe_e*n(idx_e) +1.D-9
       M30Fe = g30Fe_e*n(idx_e) +2.D-3
       M40Fe = g40Fe_e*n(idx_e) +1.5D-3
       M43Fe = g43Fe_e*n(idx_e) +3.6D-3

       !//O
       M10O = g10O_H2o*nH2o +g10O_H2p*nH2p +g10O_H*n(idx_H) +g10O_Hp*n(idx_Hj) &
            +g10O_e*n(idx_e) +8.9D-5
       M21O = g21O_H2o*nH2o +g21O_H2p*nH2p +g21O_H*n(idx_H) +g21O_Hp*n(idx_Hj) &
            +g21O_e*n(idx_e) +1.3D-10
       M20O = g20O_H2o*nH2o +g20O_H2p*nH2p +g20O_H*n(idx_H) +g20O_Hp*n(idx_Hj) &
            +g20O_e*n(idx_e) +1.8D-5
       !//Cp
       M10Cp= g10Cp_H*n(idx_H) +g10Cp_e*n(idx_e) +2.4D-6
       !//Op
       M10Op= g10Op_e*n(idx_e) +5.1D-5
       M20Op= g20Op_e*n(idx_e) +1.7D-4
       M21Op= g21Op_e*n(idx_e) +1.3D-7
       !//Sip
       M10Sip= g10Sip_e*n(idx_e) +g10Sip_H*n(idx_H) +2.2D-4
       !//Fep
       M10Fep= g10Fep_e*n(idx_e) +g10Fep_H*n(idx_H) +2.13D-3
       M21Fep= g21Fep_e*n(idx_e) +g21Fep_H*n(idx_H) +1.57D-3
       M32Fep= g32Fep_e*n(idx_e) +g32Fep_H*n(idx_H) +7.18D-4
       M43Fep= g43Fep_e*n(idx_e) +g43Fep_H*n(idx_H) +1.88D-4
       M20Fep= g20Fep_e*n(idx_e) +g20Fep_H*n(idx_H) +1.50D-9
       M30Fep= g30Fep_e*n(idx_e) +g30Fep_H*n(idx_H) +1.12D-3
       M40Fep= g40Fep_e*n(idx_e) +g40Fep_H*n(idx_H) +1.12D-3
       M31Fep= g31Fep_e*n(idx_e) +g31Fep_H*n(idx_H) +1.12D-3
       M41Fep= g41Fep_e*n(idx_e) +g41Fep_H*n(idx_H) +1.12D-3
       M42Fep= g42Fep_e*n(idx_e) +g42Fep_H*n(idx_H) +1.12D-3


       !(4)-----------------------------------------------------------------
       !//matrici di coefficienti
       AC(1,:) = (/1.d0, 1.d0, 1.d0/) 
       AC(2,:) = (/M01C+M02C, -M10C, -M20C/)
       AC(3,:) = (/-M01C, M10C+M12C, -M21C/)
       BC(:) = (/n(idx_C), 0.d0, 0.d0/)   

       ASi(1,:) = (/1.d0, 1.d0, 1.d0/)
       ASi(2,:) = (/M01Si+M02Si,-M10Si,-M20Si/)
       ASi(3,:) = (/-M01Si,M10Si+M12Si,-M21Si/)
       BSi(:) = (/n(idx_Si), 0.d0, 0.d0/)

       AFe(1,:)=(/1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
       AFe(2,:)=(/M01Fe+M02Fe+M03Fe+M04Fe,-M10Fe,-M20Fe,-M30Fe,-M40Fe/)
       AFe(3,:)=(/-M01Fe,M10Fe+M20Fe,-M21Fe,0.d0,0.d0/)
       AFe(4,:)=(/-M02Fe, -M12Fe, M20Fe+M21Fe,0.d0,0.d0/)
       AFe(5,:)=(/-M03Fe,0.d0,0.d0,M30Fe+M34Fe,-M43Fe/)
       BFe(:) = (/n(idx_Fe), 0.d0, 0.d0, 0.d0, 0.d0/)

       AO(1,:) = (/1.d0, 1.d0, 1.d0/) 
       AO(2,:) = (/M01O+M02O,-M10O,-M20O/)
       AO(3,:) = (/-M01O,M10O+M12O,-M21O/)
       BO(:) = (/n(idx_O), 0.d0, 0.d0/)

       ACp(1,:) = (/1.d0, 1.d0/)
       ACp(2,:) = (/M01Cp,-M10Cp/)
       BCp(:) = (/n(idx_Cj), 0.d0/)

       AOp(1,:)=(/1.d0, 1.d0, 1.d0/)
       AOp(2,:)=(/M01Op+M02Op,-M10Op,-M20Op/)
       AOp(3,:)=(/-M01Op,M10Op+M12Op,-M21Op/)
       BOp(:)=(/n(idx_Oj),0.d0,0.d0/)

       ASip(1,:)=(/1.d0, 1.d0/)
       ASip(2,:)=(/M01Sip,-M10Sip/)
       BSip(:)=(/n(idx_Sij), 0.d0/)

       AFep(1,:)=(/1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
       AFep(2,:)=(/M01Fep+M02Fep+M03Fep+M04Fep,-M10Fep,-M20Fep,-M30Fep,-M40Fep/)
       AFep(3,:)=(/-M01Fep,M10Fep+M20Fep,-M21Fep, 0.d0, 0.d0/)
       AFep(4,:)=(/-M02Fep, -M12Fep, M20Fep+M21Fep,0.d0,0.d0/)
       AFep(5,:)=(/-M03Fep,0.d0,0.d0,M30Fep+M34Fep,-M43Fep/)
       BFep(:) = (/n(idx_Fej),0.d0,0.d0,0.d0,0.d0/)


       !(5)------------------------------------------------------------------
       call mydgesv(AC,BC)
       call mydgesv(ASi,BSi)
       call mydgesv(AFe,BFe)
       call mydgesv(AO,BO)
       call mydgesv(ACp,BCp)
       call mydgesv(ASip,BSip)
       call mydgesv(AFep,BFep)
       call mydgesv(AOp,BOp)

       !(6)------------------------------------------------------------------
       cool = cool + (BC(2)*7.9D-8*24. +BC(3)*2.1D-14*63. +BC(3)*2.7D-7*39.)*kb
       cool = cool + (BSi(2)*8.4D-6*110. +BSi(3)*2.4D-10*320. &
            +BSi(3)*4.2D-5*210.)*kb
       cool = cool + (BFe(2)*2.5D-3*600.+BFe(3)*(1.D-9*1000.+1.6D-3*420.)&
            +BFe(4)*2.D-3*9900. +BFe(5)*(1.5D-3*11000.+3.6D-3*640.))*kb
       cool = cool + (BO(2)*8.5D-5*230.+BO(3)*(1.3D-10*330.+1.8D-5*98.))*kb
       cool = cool + (BCp(2)*2.4D-6*1.259D-14)
       cool = cool + (BOp(2)*5.1D-5*3.9D4 +BOp(3)*(1.7D-4*3.9D4+1.3D-7*30.))*kb
       cool = cool + (BSip(2)*2.2D-4*410.)*kb
       cool = cool + (BFep(2)*2.13D-3*553.58 +BFep(3)*(1.5D-9*960.59&
            +1.57D-3*407.01) + &
            BFep(4)*(1.12D-3*1241.16 +1.12D-3*687.58 +1.12D-3*280.57) + &
            BFep(5)*(1.12D-3*1405.76 +1.12D-3*852.18 +1.12D-3*445.17 &
            +1.88D-4*164.6))*kb
    else
       cool = 0.d0
    end if


    !check errore
    if(cool<0.d0) then
       print *,"ERROR!!! metal cooling <0",cool
       do i=1,size(n)
          print '(I3,E11.3)',i,n(i)
       end do
       stop
    end if
    cooling_Z = cool
  end function cooling_Z

  !*********************************
  subroutine mydgesv(A,B)
    real*8::A(:,:),B(:)
    integer::n,info,i
    integer,allocatable::ipiv(:),tmp(:)
    n=size(B)
    allocate(tmp(n))
    allocate(ipiv(n))
    call dgesv(n,1,A,n,ipiv,B,n,info)
    if(info .ne. 0) then
       print *,"ERRORE!!! problemi con dgesv, info: ",info
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

#IFKROME_useHeatingPhoto
  !**************************
  function photo_heating(n)
    use krome_commons
    real*8::photo_heating,n(:),n0
    n0 = 1d-3 !density for fake opacity (see KROME paper)
    photo_heating = 0.d0
#KROME_photo_heating
  end function photo_heating
#ENDIFKROME

#IFKROME_useHeatingA
  !H2 FORMATION HEATING
  !OMUKAI 2000
  !UNITS=erg/cm3/s
  !Following Hollenbach & McKee (1979), we assume the heat deposited 
  !per a formed H2 is: 
  !3.53*(1+ncr/n)^−1 eV for H2 formation by H− process (k_HM) 
  !1.83*(1+ncr/n)^−1 eV for H2 formation by H2+ process (k_H2j)
  !4.48*(1+ncr/n)^−1 eV for H2 formation by the three-body 
  !reactions (k_3B, k_3B1)
  !NOTE: THE LAST TERM IS A COOLING SOURCE TERM FROM k_3B2 COSA VUOI FARE?
  !*******************************
  function heatingA(n, Tgas)
    use krome_constants
    use krome_commons
    implicit none
    real*8::heatingA, n(:), Tgas
    real*8::k_3B2, k_HM, k_H2j, k_3B, k_3B1
    real*8::h2heatfac,H2delta,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2
    real*8::Te,invTe,dd

    dd = sum(n(1:nmols))

    Te = Tgas*boltzmann_eV
    invTe = 1.d0/te

    heatingA = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**(2.d0))
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    k_3B  = 5.5d-29/Tgas  !3H --> H2 + H
    k_3B1 = k_3B/8.0d0    !2H + H2 --> H2 + H 
    k_3B2 = (6.5d-7/Tgas**0.5d0)*exp(-5.2d4/Tgas)&
         *(1.d0-exp(-6.0d3/Tgas)) !H2 + H --> 3H

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    H2delta = n(idx_H)*(4.48d0*k_3B*n(idx_H)**(2.d0) &
         + 4.48d0*(k_3B1)*n(idx_H)*n(idx_H2)*2.0 &
                                !+ 3.53d0*k_HM*n(idx_Hk) &
                                !+ 1.83*k_H2j*n(idx_H2)
         )*h2heatfac &
         - n(idx_H)*(4.48d0*k_3B2*n(idx_H2))


    if(H2delta.gt.0.0d0)then
       heatingA = H2delta*eV_to_erg  !erg/cm3/s
    else
       heatingA = 0.0d0
    endif

  end function heatingA
#ENDIFKROME

#IFKROME_useHeatingCompress
  !***********************
  function heat_compress(n, Tgas)
    use krome_user_commons
    use krome_commons
    use krome_constants
    real*8::heat_compress,n(:), eint, gamma
    real*8::dd, p2d, Tgas
    gamma = 5./3.
    dd = sum(n(1:nmols))
    eint = (1.0/(gamma-1.0d0))*(boltzmann_erg*Tgas)/(1.22*p_mass) !erg/g 
    p2d = (gamma-1.0d0)*dd*eint  !erg/g/cm3

    !COMPRESSIONAL HEATING
    heat_compress = (p2d)*(p_mass/user_tff) !erg/s/cm3
  end function heat_compress
#ENDIFKROME

  !************************************
  subroutine plot_cool(n)
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_CEN,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP,cool_Z
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
       cool_CEN = 0.d0
       cool_Z = 0.d0
#IFKROME_useCoolingH2
       cool_H2 = cooling_H2(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingH2GP
       cool_H2GP = cooling_H2GP(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingCEN
       cool_CEN = cooling_CEN(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingHD
       cool_HD = cooling_HD(n(:),Tgas)
#ENDIFKROME
#IFKROME_useCoolingZ
       cool_Z = cooling_Z(n(:),Tgas)
#ENDIFKROME
       cool_tot = cool_H2 + cool_CEN + cool_HD + cool_Z
       cool_totGP = cool_H2GP + cool_CEN + cool_HD + cool_Z
       write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, cool_CEN, &
            cool_HD, cool_H2GP, cool_Z
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool


end module KROME_cooling
