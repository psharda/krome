module krome_subs
contains

#KROME_header

  !************************
  !compute reaction rates
  function coe(n)
    use krome_commons
    use krome_user_commons
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,t,t3,n(nspec)
    real*8::logT,lnT,Te,lnTe,T32,invT,invTe,sqrTgas,invsqrT32,sqrT32
    real*8::Tgas2,Tgas3,Tgas4,T0,T02,T03,T04,T0inv,T4
    integer::i
    !Tgas is in K
    Tgas = max(n(idx_Tgas), 2.73d0)

#KROME_Tshortcuts

    k(:) = 1.d-40 !inizialize coefficients

#KROME_krates

    !set a minumum value for the rates
    !do i=1,nrea
    !   k(i) = max(k(i),1d-40) 
    !end do

    coe(:) = k(:)!set coefficients to return variable
  end function coe

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

#KROME_masses

  end function get_mass

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

#KROME_names

  end function get_names

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
       !when found store and break loop
       if(trim(names(i))== trim(name)) then
          get_index = i !store index
          exit
       end if
    end do
    
    !error if species not found
    if(get_index<0) then
       print *,"ERROR: can't find the index of ",name
       stop
    end if

  end function get_index

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

#KROME_charges

  end function get_charges
  
  !************************
  !get species charges
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

#KROME_reaction_names

  end function get_rnames

  !*****************************
  function revKc(Tgas,ridx,pidx)
    real*8::revKc,Tgas
    integer::ridx(:),pidx(:),i

    revKc = 0.d0

    do i=1,size(pidx)
       revKc = revKc + revHS(Tgas,pidx(i))
    end do

    do i=1,size(ridx)
       revKc = revKc - revHS(Tgas,ridx(i))
    end do

  end function revKc
  
  !*****************************
  function revHS(Tgas,idx)
    use krome_commons
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
#KROME_var_reverse
    integer::idx

    p1(:,:) = 0.d0
    p2(:,:) = 0.d0
    Tlim(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)

#KROME_kc_reverse

    if(Tlim(idx,2)==0.d0) then
       revHS = 0.d0
       return
    end if

    if(Tlim(idx,1).le.Tgas .and. Tgas.le.Tlim(idx,2)) p(:) = p1(idx,:)
    if(Tlim(idx,2)<Tgas .and. Tgas.le.Tlim(idx,3)) p(:) = p2(idx,:)
    H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
         p(5)*Tgas4*0.2d0 + p(6)*invT
    S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
         p(5)*Tgas4*0.25d0 + p(7)

    revHS = H - S

  end function revHS

  !*****************************
  subroutine load_arrays()
    use krome_commons

#KROME_implicit_arrays

  end subroutine load_arrays

end module krome_subs
