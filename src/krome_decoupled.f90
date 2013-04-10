module krome_decoupled
contains

  !**********************************************************
  ! EVALUATE THE ENERGY VARIATION DUE TO COOLING AND HEATING
  ! UNITS: erg/grams
  ! ABUNDANCES for heat/cool ARE NEEDED IN NUMBER DENSITY
  subroutine krome_evaluate_edot(x,temp,rhogas,edot)
    use krome_commons
    use krome_subs
    use krome_cooling
    use krome_heating
    implicit none
    integer::i
    real*8::temp,edot,rhogas
    real*8::dt,x(nmols),mass(nspec),n(nspec)

    !edot = Lambda - Gamma (in units of erg/g)
    mass(:) = get_mass() !get masses
    !compute densities from fractions
    do i = 1,nspec
       if(mass(i)>0.d0) n(i) = rhogas * x(i) / mass(i)
    end do

    edot = heating(n(:), temp) - cooling(n(:), temp)

  end subroutine krome_evaluate_edot
  !*********************************************************

  !*********************************************************
  ! EVALUATE THE TIMESTEP BASED ON THE ENERGY VARIATION
  ! KEEPING THE CHANGE BELOW 10%
  subroutine krome_timestep(energy,tg,edot,dt,tot,dt_cl)
    implicit none
    real*8::energy,edot,tg
    real*8::dt,dt_cl,tot,olddtit

    ! energy and edot in erg/grams
    if(tg.lt.10.0d0 .and. edot.lt.0.0d0) edot = 1.d-30
    dt_cl = min(real(abs(0.1d0*energy/edot)),&
         dt-tot, dt)
  end subroutine krome_timestep
  !*********************************************************

  !*********************************************************
  ! ESTABLISH THE EXIT STATUS FOR THE SUBCYCLES
  subroutine krome_timestep_exit(nx,ttot,dt,dt_cool,ttmin)
    implicit none
    integer::nx,i
    real*8::dt,ttmin
    real*8::ttot(nx),dt_cool(nx)

    ttmin = 1.d30                                                      
    do i = 1, nx
       ttot(i) = min(ttot(i) + dt_cool(i), dt)
       if (ttot(i).lt.ttmin) ttmin = ttot(i)
    enddo

  end subroutine krome_timestep_exit
  !**********************************************************

end module krome_decoupled
