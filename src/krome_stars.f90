module krome_stars
#IFKROME_useStars
  use krome_commons
  implicit none
  integer,parameter::nzp=26
  real*8::pow13(nzp),pow23(nzp),pow43(nzp)
  real*8::pow53(nzp),pow186(nzp),pow158(nzp)
contains


  !********************
  !initialize the common arrays for stellar physics
  subroutine stars_init()
    use krome_subs
    implicit none
    integer::i
    do i=1,nzp
       pow13(i) = (1.0d0*i)**(1./3.)  !! ^1/3
       pow23(i) = (1.0d0*i)**(2.0/3.0) !! ^2/3
       pow43(i) = (1.0d0*i)**(4.0/3.0) !! ^4/3
       pow53(i) = (1.0d0*i)**(5.0/3.0) !! ^5/3
       pow186(i) = (1.0d0*i)**(1.860)  !! ^1.860
       pow158(i) = (1.0d0*i)**(1.580)  !! ^1.580
    end do

  end subroutine stars_init


  !*****************************
  !adapt rate coefficients for nuclear reactions in stars,
  ! using a template screening function. more details in 
  ! stars_screen function
  function stars_coe(n,rho,Tgas,y,zz)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8::n(:),Tgas,k(nrea),stars_coe(nrea),rho,ko(nrea)
    real*8::scr12,scr23,y(:)
    integer::zz(:),z1,z2,i,z12,zo(nspec)

    zo(:) = get_zatoms()
    ko(:) = coe(n(:))

    !scale reactions using screening
    do i=1,nrea
       z1 = zo(arr_r1(i)) !first reactant atomic number
       z2 = zo(arr_r2(i)) !second reactant atomic number
       k(i) = ko(i) &
            * stars_screen(Tgas,rho,y(:),zz(:),z1,z2)
    end do

#KROME_stars_3body

    stars_coe(:) = k(:)

  end function stars_coe

  !*************************
  !screening function. WARNING: the screening provided
  ! is provided as a template since is suited only for 
  ! ionized medium, two-body interaction, and nuclear-only 
  ! interaction.
  function stars_screen(Tgas,rho, n, zz, z1, z2)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::stars_screen
    real*8::Tgas,rho,n(:),ZY,T13,Atmp,pmol,lamb0,lamb12
    real*8::lamb0i,Zt058,Zm028,Z158m,H12i,lamb0s,Zm13
    real*8::a,b,c,d,CHK,H12,fer,theta,Zt,Zm,f(size(n))
    real*8::H12s,zz158(size(n))
    integer::z1,z2,i,zz(:)
    real*8,parameter::H12max=2d2

    !screen following
    ! GRABOSKE: Graboske, DeWitt, Grossman, and Cooper, 1973, ApJ, 181, 457
    !           DeWitt, Graboske, and Cooper, 1973, ApJ, 181, 439
    ! ITOH:     Itoh, Totsuji, and Ichimaru, 1977, ApJ, 218, 477
    !           Itoh, Totsuji, Ichimaru, and DeWitt, 1979, ApJ, 234, 1079

    zz158(:) = zz(:)**1.580

    !ITOH
    ZY = sum(zz(:)*n(:))
    Atmp = 1.442250 / (7.567682d24 * rho * ZY)**(1./3.) 
    T13 = Tgas**(1./3.)

    !strong screening
    a   = pow13(z1) * Atmp 
    b   = pow13(z2) * Atmp 
    c  = (z1 * z2 * 2.307121d-19) / (0.5 * (a + b) * boltzmann_erg * Tgas) 
    d = z1 * z2 * 5351.69 / (pow13(z1 + z2) * T13)
    CHK  = 3.0 * c / d

    if(CHK >= 0.20d0) then 
       H12 = 1.25 * c - 0.095 * d * CHK**2
       stars_screen = exp(H12) 
       return 
    end if

    !GRABOSKE & DeWitt
    !Theta
    fer = 5.4885d7 * rho * (1.0d0 + n(idx_H)) /(Tgas * sqrt(Tgas))
    if(fer <= 3.0d0) then
       theta = 1.0d0 / (1.0d0 + 0.39716 * fer - 0.00929 * (fer**2) )
    else
       theta = 1.1447d0 / ( (fer**(2./3.)) * (1.0d0 + 0.28077 / fer) )
    end if

    pmol = 1.d0/sum(n(:))
    Zm = sum(ZZ(:) * n(:))
    Zt = sum(ZZ(:) * n(:) * (ZZ(:) + theta))

    f(:) = n(:) * pmol
    Zm = Zm * pmol
    Zt = sqrt(Zt * pmol)

    lamb0 = 1.88d8 * sqrt(rho / (Tgas * Tgas * Tgas * pmol))
    lamb12 = lamb0 * Z1 * Z2 * Zt

    if(lamb12.le.1d-1) then
       !weak screening
       stars_screen = exp(lamb12)
    elseif(lamb12.le.5.d0) then
       !intermediate screening
       lamb0i = lamb0**(0.860)
       a = pow186(Z1 + Z2) - pow186(Z1) - pow186(Z2)
       Zt058 = Zt**(0.580)
       Zm028 = Zm**(0.280)

       Z158m = sum(f(:) * ZZ158(:))

       H12i = 0.38d0 * a * lamb0i * Z158m / (Zt058 * Zm028)

       if(lamb12.le.2.d0) then
          !intermediate GRABOSKE
          stars_screen = exp(H12i)
       else
          !intermediate-strong GRABOSKE
          !minimum between intermediate and strong

          !! schermaggio forte !!
          lamb0s = lamb0**(2.0/3.0)
          a = pow53(Z1 + Z2) - pow53(Z1) - pow53(Z2)
          b = pow43(Z1 + Z2) - pow43(Z1) - pow43(Z2)
          c = (pow23(Z1 + Z2) - pow23(Z1) - pow23(Z2) ) / lamb0s
          Zm13 = Zm**(1.0/3.0)
          H12s = 0.624d0 * Zm13 * &
               (a + 0.316d0 * b * Zm13 + 0.737d0 * c / Zm13 ) * lamb0s
          !controls exceeding
          if(H12s > H12max) H12s = H12max
          H12 = min(H12i, H12s)
          stars_screen = exp(H12)
       end if
    else
       ! GRABOSKE !!
       ! strong screening
       lamb0s = lamb0**(2.0/3.0)
       a = pow53(Z1 + Z2) - pow53(Z1) - pow53(Z2)
       b = pow43(Z1 + Z2) - pow43(Z1) - pow43(Z2)
       c = (pow23(Z1 + Z2) - pow23(Z1) - pow23(Z2) ) / lamb0s
       Zm13 = Zm**(1.0/3.0)
       H12 = 0.624d0 * Zm13 * (a + 0.316d0 * b * Zm13 + 0.737d0 * c / Zm13 ) * &
            & lamb0s

       !controls exceeding
       if(H12 > H12max) H12 = H12max
       stars_screen = exp(H12)
    end if

  end function stars_screen


  !*****************************
  !computes energies for nuclear reactions using Qeff and flux
  function stars_energy(n,rho,Tgas,k)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8::stars_energy(nrea),qeff(nrea),flux(nrea)
    real*8::k(nrea),n(:),rho,Tgas
    qeff(:) = get_qeff()
    
#KROME_stars_energy

    stars_energy(:) = flux(:) * qeff(:)

  end function stars_energy

#ENDIFKROME

end module krome_stars
