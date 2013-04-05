!****************SEDOV-LIKE SHOCK 1D LAGRANGIAN TEST***************
! REFER TO: Bodenheimer P., Laughlin G., Rozyczka M., Yorke H.    
!           2006, Numerical Methods in Astrophysics:         
!           An Introduction. Series in Astronomy and Astrophysics, 
!                             Taylor&Francis 
!*******************************************************************

!The scheme of the following f90 file is:
!module commons
!module subs
!program sedov (main)

!***************START MODULES NEEDED BY THE TEST-CODE***************
!*****COMMONS*****!
module commons
  use krome_commons
  integer,parameter::nx=100, il=20, nsp=nmols
  real*8,parameter:: pi43=4.d0/3.d0*3.14, spy=365.*24.*3600.d0
  real*8::cflfactor,dfactor,gamma,dt12,dt,t, fmratio, efac, tmax, Tgas(nx)
  real*8::radius, rhogas, mass0, esed

  real*8::dm12(nx), rho(nx), eps(nx), p(nx), u(nx), fm(nx), r(nx), v(nx)
  real*8::at12(nx), rold(nx), a(nx), dr12(nx), dm(nx), w(nx), ak12(nx)
  real*8::x(nx,nsp),massr(nx)
end module commons

!******SUBS***********! 
module subs
contains
  !*****INITIALIZATION OF THE CHEMICAL ABUNDANCES*******!
  subroutine initchem()
    use commons
    use krome_user

    integer::ix
    real*8:: xx(nsp)
    x(:,:) = 0.d0

    !   INITIALIZATION IN MASS FRACTION
    do ix = 1, nx
       x(ix,KROME_idx_H)     = 0.9225d0  !H
       x(ix,KROME_idx_E)     = 0.d0 !1.0d-4    !E
       x(ix,KROME_idx_Hj)    = 0.d0 !1.0d-4    !H+
       x(ix,KROME_idx_D)     = 1.0d-20   !D
       x(ix,KROME_idx_Dj)    = 0.d0 !1.0d-20   !D+
       x(ix,KROME_idx_HE)    = 0.0972d0  !He
       x(ix,KROME_idx_HEj)   = 0.d0 !1.0d-20   !He+
       x(ix,KROME_idx_H2j)   = 0.d0 !1.0d-20   !H2+
       x(ix,KROME_idx_H2)    = 1.0d-5    !H2
       x(ix,KROME_idx_HD)    = 1.0d-8    !HD
       x(ix,KROME_idx_Hk)    = 0.d0 !1.0d-20   !H-
       x(ix,KROME_idx_HEjj)  = 0.d0 !1.0d-20   !He++

       x(ix,:) = x(ix,:) / sum(x(ix,:))
    end do

    do ix = 1, nx
       xx(:) = x(ix,:) 
       !write(32,'(99E11.3)') xx(:)
    end do

  end subroutine initchem

  !******TIME STEP EVALUATION********!
  subroutine tmsc
    use commons
    real*8::dtc,dtd
    integer::ix

    dtc = 1.d99
    do ix = 1, nx-1
       dtc = min(dtc, dr12(ix) / (abs(u(ix)) + sqrt(gamma*eps(ix))))
    end do

    !courant
    dtc  = cflfactor * dtc
    if (t+dtc.gt.tmax) dtc = tmax-t 

    !diffusion limit
    dtd = 1.d-99
    do ix = 1, nx-1
       dtd = max(dtd, abs(at12(ix+1)*u(ix+1)-at12(ix)*u(ix)) &
            /(v(ix+1) - v(ix)))
    end do

    dtd = .5d0/dtd/dfactor
    dtc = min(dtc,dtd)
    dt   = .5d0 * (dt12 + dtc)
    dt12 = dtc
  end subroutine tmsc

  !*****HYDRODINAMIC INITIAL CONDITIONS********!
  subroutine inicond
    use commons
    use krome_constants
    implicit none
    real*8::dmej,rhoej,rhoamb,dmamb,rSN
    integer::ix

    !Courant factor, must be smaller than 1
    cflfactor = .1d0

    fmratio = .1d0

    !diffusion factor, must be greater that 1
    dfactor = 1.5d0
    dfactor = dfactor**2

    rSN = radius !cm

    !high pressure ejecta
    rhoej  =  ((dble(nx)/dble(il))**3 - 1.d0) * fmratio * rhogas
    dmej   =  pi43 * (dble(il)/dble(nx) * rSN)**3 * rhoej / dble(il)
    efac = 1d6 * boltzmann_erg / p_mass

    print *,"***********************"
    print *,"rhoej (g/cm3):",rhoej
    print *,"dmej (g):",dmej
    print *,"epsej (erg/g):",efac

    do ix = 1, il
       dm12(ix) = dmej
       rho (ix) = rhoej
       eps (ix) = efac
       p   (ix) = (gamma - 1d0) * rho(ix) * eps(ix)
       u   (ix) = 0.d0
    end do

    print *,"***********************"
    !low pressure ambient medium
    rhoamb = rhogas
    dmamb  = pi43 * rSN**3 * (1.d0-(dble(il)/dble(nx))**3) / dble(nx-il) &
         * rhoamb
    efac = 1d1 * boltzmann_erg / p_mass
    print *,"rhoamb (g/cm3):",rhoamb
    print *,"dmamb (g):",dmamb
    print *,"epsamb (erg/g):",efac

    do ix = il + 1, nx
       dm12(ix) = dmamb
       rho (ix) = rhoamb
       eps (ix) = efac
       p   (ix) = (gamma - 1.d0) * rho(ix) * eps(ix)
       u   (ix) = 0.d0
    end do

    fm(1) = dm12(1) 
    do ix = 2, nx
       fm (ix) = fm(ix-1) + dm12(ix)
    end do
    do ix = 2, nx
       dm(ix) = 0.5d0 * (dm12(ix) + dm12(ix-1))
    end do


    r(1) = 0.d0
    v(1) = 0.d0
    do ix = 2, nx
       v(ix) = v(ix-1) + dm12(ix-1) / rho(ix-1)
       r(ix) = (v(ix)/pi43)**(1./3.)
       a(ix) = 4.d0 * pi * r(ix)**2
    end do
    do ix = 1, nx-1
       dr12(ix) = r(ix+1) - r(ix)
    end do

    print *,"***********************"
    print *,"Rtot:",r(nx)
    print *,"Vtot:",v(nx)
    print *,"Mtot:",sum(dm(:))

    do ix = 1,nx
       write(22,'(99E11.3)') r(ix), dr12(ix), dm12(ix), eps(ix), rho(ix)
    end do
    t = 0.d0
    print *,"***********************"

  end subroutine inicond
end module subs
!*********************END MODULES********************************!


!*****************START PROGRAM SEDOV****************************!
program sedov

  use krome_main
  use krome_user
  use krome_constants
  use krome_subs
  use krome_photo
  use commons
  use subs
  implicit none
  integer::nsteps, istep, ix, ndump,i
  real*8::q,aux,ethe,ekin,etot,etot0
  real*8::umax, rhomax, pmax, emax, wmax, epsi
  real*8::xx(nsp)

  !     This Lagrangian code follows the adiabatic expansion 
  !     of a hot spherical bubble in a uniform ambient medium.
  !     The present setup is ajusted to spherical geometry.
  !     It includes primordial chemistry evolution
  !     and H2/HD/CEN cooling

  !open the file in which the final model will be stored 
  open(16,file='lh1.out',form='formatted',status='replace')

  !open the file in which values of thermal, kinetic and 
  ! energy will be stored at 10th time-step
  open(17,file='lh2.out',form='formatted',status='replace')

  call krome_init()
  call krome_init_photo()

  gamma = 5.d0/3.d0
  nsteps = int(1e4) !maximum number of time-steps
  ndump = 10 !dump interval
  Tgas(:) = 1d3 !default gas temp (K)  
  tmax = spy * 5d3 !maximum integration time (s)
  q = 2d0 !artificial viscosity parameter
  radius = 3.08568025d18 * 1d0 !radius of the box (cm)
  rhogas = 1d-24 !gas density (g/cm3)

  xx(:) = 0.d0

  !write input to check
  write(*,221) "nsteps:", nsteps !number of time steps
  write(*,222) "tmax (s):", tmax   !maximum time
  write(*,222) "Tgas (K):", sum(Tgas)/dble(nx)   !gas Temp
  write(*,222) "q:", q      !artificial viscosity parameter
  write(*,222) "radius (cm):", radius !radius of the object
  write(*,222) "density (g/cm3):", rhogas !density of the object
  write(*,222) "density (1/cm3):", rhogas/p_mass !number density of the object

221 format(a20,I11)
222 format(a20,E11.3)


  call inicond()  !set initial conditions
  call initchem() !init chemistry, abundances in mass fraction

  !begin loop over time steps
  do istep = 1, nsteps

     !determine time step
     call tmsc

     !update velocities (eqn 6.42)
     do ix = 2, nx
        u(ix) = u(ix) - a(ix) * (p(ix)-p(ix-1)) * dt / dm(ix) &
             - .5d0 * (w(ix) * (3.d0*ak12(ix)-a(ix)) &
             - w(ix-1) * (3.d0*ak12(ix-1)-a(ix))) &
             * dt/dm(ix) 
        u(ix) = u(ix) - gravity * massr(ix) / r(ix)**2 * dt
     end do

     u(1) = 0.d0

     !update radii, surfaces and volumes (eqn 6.41_2)
     rold(:) = r(:)
     r(:) = rold(:) + u(:) * dt12

     if(minval(r)<0.d0) then
        print *,"r(ix)<0"
        stop
     end if

     do ix = 1, nx-1
        dr12(ix) = r(ix+1) - r(ix)
     end do

     at12(:) = 4.d0 * pi * (.5d0 * (r(:)+rold(:)))**2
     a(:) = 4.d0 * pi * r(:)**2 
     v(:) = pi43 * r(:)**3

     do ix= 1, nx-1
        ak12(ix) = .5d0 * (at12(ix+1) + at12(ix))
     end do

     !update densities (6.41_3)
     do ix = 1, nx-1
        rho(ix)  = dm12(ix) / (v(ix+1)-v(ix))
     end do
     rho(nx)  = rho(nx-1)

     !update shells integral mass 
     massr(1) = 0.d0
     do ix = 2, nx
        massr(ix) = massr(ix-1) + pi43 * (r(ix)**3 -r(ix-1)**3) * rho(ix)
     end do

     !artificial viscosity
     do ix = 1, nx-1
        w(ix)  = -q**2 * rho(ix) * abs(u(ix+1)-u(ix)) &
             * (u(ix+1) * (1.d0 - at12(ix+1)/3.d0/ak12(ix)) &
             - u(ix) * (1.d0 - at12(ix)/3.d0/ak12(ix))) 
     end do

     do ix = 1, nx-1
        if(u(ix+1) > u(ix)) w(ix) = 0.d0 !no viscosity for postive velocity divergence
     end do

     !update internal energies and pressures (6.43)
     do ix = 1, nx-1
        aux = eps(ix) - p(ix) &
             * (at12(ix+1) * u(ix+1) &
             - at12(ix) * u(ix)) * dt12 / dm12(ix)
        p(ix) = .5d0 * (p(ix) + (gamma-1.d0) * rho(ix) * aux)
     end do

     do ix = 1, nx-1
        eps(ix) = eps(ix) - p(ix) &
             * (at12(ix+1) * u(ix+1) &
             - at12(ix) * u(ix)) * dt12/dm12(ix)
     end do

     !contribution from artificial viscosity
     do ix = 1, nx-1
        eps(ix) = eps(ix) - .5d0 * w(ix) * dt12/dm12(ix) &
             * (u(ix+1)*(3.d0 * ak12(ix) - at12(ix+1)) &
             - u(ix)*(3.d0 * ak12(ix) - at12(ix)))
        Tgas(ix) = (gamma - 1.d0) * eps(ix) / boltzmann_erg * p_mass
     end do

     !DO CHEMISTRY: CALL KROME PACKAGE
     !REMEMBER xx in fraction, density in g/cm-3
     do ix = 1, nx
        xx(:) = x(ix,:)
        rhogas = rho(ix)
        call krome(xx(:),rhogas,Tgas(ix),dt)
        x(ix,:) = xx(:)
     enddo
     !update internal energy with new temperature
     eps(:) = Tgas(:) / (gamma - 1.d0) / p_mass * boltzmann_erg
     !*********END******************************!

     !update pressure
     do ix = 1, nx-1
        p(ix) = (gamma-1.d0) * rho(ix) * eps(ix)
     end do
     !update last shell
     p(nx) = p(nx-1)
     eps(nx) = eps(nx-1)

     !end time step

     t = t + dt12  

     !check energy conservation
     ethe = 0.d0
     ekin = 0.d0
     do ix = 2, nx-1
        ethe = ethe + (eps(ix)) * dm(ix)
        ekin = ekin + .5d0 * (.5d0 * (u(ix+1)+u(ix)))**2 * dm(ix)
     end do
     etot = ethe + ekin
     if (istep == 1) etot0 = etot
     etot = etot / etot0
     ethe = ethe / etot0
     ekin = ekin / etot0

     if (mod(istep,10).eq.0)  write(17,200) istep, t, etot, ethe, ekin
     if (mod(istep,ndump).eq.0) write(* ,100) istep, t/spy, etot, ethe, ekin

     if (mod(istep,ndump).eq.0) then
        do ix = 1, nx
           write(23,'(I8,99E11.3)') ix, t, r(ix), rho(ix), abs(u(ix)), &
                p(ix), Tgas(ix)
        end do
        write(23,*)
     end if

     !end loop over time steps
     if (t.ge.tmax) exit
  end do

  !dump final conditions
  !dump chemistry
  do ix = 1, nx
     !get number density from mass fraction
     xx(1:nsp) = x(ix,:)!*rho(ix)/spec_mass(:)
     write(24,'(I8,99E11.3)') ix, t, r(ix), rho(ix), abs(u(ix)), p(ix), Tgas(ix)
     write(34,'(I8,999E11.3)') ix, t, r(ix), xx(1:nsp)
  end do
  write(24,*)


  !final printout

  umax = 0.d0
  rhomax = 0.d0
  pmax = 0.d0
  emax = 0.d0
  wmax = 0.d0

  do ix = 1, nx
     if (u(ix) > 0) umax  = max(u(ix),umax)
     rhomax = max(rho(ix),rhomax)
     pmax   = max(p(ix),pmax)
     emax   = max(eps(ix),emax)
     wmax   = max(w(ix),wmax)
  end do

  epsi = 1.e-20

  do ix=1,nx
     u(ix)   = u(ix) / (epsi+umax)
     rho(ix) = rho(ix) / (epsi+rhomax)
     p(ix)   = p(ix) / (epsi+pmax)
     eps(ix) = eps(ix) / (epsi+emax)
     w(ix)   = w(ix) / (epsi+wmax)
  end do

  do ix = 1, nx
     write(16,101) ix,fm(ix),r(ix)/r(nx),u(ix),rho(ix),p(ix),eps(ix),w(ix)
  end do

100 format(1x,i8,'; t:',1pe10.2,';  etot:',1pe10.2,';  eth:',1pe10.2,';  ekin:',1pe10.2)
101 format(1x,i8,1p7e12.4) 
200 format(1x,i8,1p5e11.3)
  close(16)
  close(17)


  write(*,*) "finish!"
end program sedov
