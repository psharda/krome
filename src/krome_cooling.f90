module KROME_cooling
contains

#KROME_header


#IFKROME_useCooling
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
  end function cooling
#ENDIFKROME

  !*******************************
  function heating(n, Tgas)
    implicit none
    real*8::n(:), Tgas
    real*8::heating 
    !total heating erg/cm3/s
    heating = 0.d0
#IFKROME_useHeating
    heating = heatingA(n(:), Tgas) + heat_compress(n(:), Tgas)
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

#IFKROME_useHeating
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
    real*8::cool_CEN,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
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
       cool_tot = cool_H2 + cool_CEN + cool_HD
       cool_totGP = cool_H2GP + cool_CEN + cool_HD
       write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, cool_CEN, &
            cool_HD, cool_H2GP
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool


end module KROME_cooling
