module krome_reduction
contains


#IFKROME_useTopology
  subroutine krome_authub(threshold)
    use krome_commons
    use krome_subs
    implicit none

    integer::i,j,k,iter,nused,uspec(nspec),rp(3+4)
    real*8::hub(nspec),aut(nspec)
    real*8::hubx(nspec),autx(nspec)
    real*8::rms,naut,nhub,thold
    real*8,optional::threshold
    character*16::names(nspec)
    !check if commons array are init
    if(maxval(arr_p1)==0) then
       print *,"ERROR: reactants/products arrays not initialized!"
       stop
    end if

    thold = 1d0 !default value of threshold
    if(present(threshold)) thold = threshold

    !initialize names and hub/auth values
    names(:) = get_names()
    hub(:) = 1.d0 / nspec
    aut(:) = 1.d0 / nspec

    iter = 0 !number of iterations
    do 
       iter = iter + 1
       hubx(:) = 0.d0
       autx(:) = 0.d0
       !loop on reactions
       do i=1,nrea
          !reactants+products indexes
          rp = (/arr_r1(i), arr_r2(i), arr_r3(i), &
               arr_p1(i), arr_p2(i), arr_p3(i), arr_p3(i)/)
          !loop on reactants
          do j=1,3
             if(rp(j)==idx_dummy) cycle
             !loop on products
             do k=4,7
                if(rp(k)==idx_dummy) cycle
                hubx(rp(j)) = hubx(rp(j)) + aut(rp(k))
                autx(rp(k)) = autx(rp(k)) + hub(rp(j))
             end do
          end do
       end do
       !normalize
       hubx(:) = hubx(:) / sum(hubx)
       autx(:) = autx(:) / sum(autx)
       !compute rms
       rms = sqrt(sum((hubx(:)-hub(:))**2)/nspec)
       rms = 0.5d0 * (rms + sqrt(sum((autx(:)-aut(:))**2)/nspec))
       !store old values
       hub(:) = hubx(:)
       aut(:) = autx(:)
       if(rms<1d-10.or.iter>5000) exit
    end do

    print *,"Authub convergence found after iterations:",iter
    print *,"with RMS:",rms
    print *,""

    !list and select usable species
    print *,"Listing species"
    print '(a5,a10,2a12)', "idx", "name", "auth", "hub"
    nused = 0
    uspec(:) = 0

    do i=1,nspec
       naut = aut(i)*1d2/maxval(aut)
       nhub =  hub(i)*1d2/maxval(hub)
       if(naut>thold .or. nhub>thold) then
          nused = nused + 1
          uspec(i) = 1
          print '(I5,a10,2F12.4,a4)', i," "//names(i), naut, nhub,""
       else
          print '(I5,a10,2F12.4,a4)', i," "//names(i), naut, nhub," XX"
       end if
    end do
    print *,"used species:",nused
    print *,"XX above threshold",thold
    print *,""


    arr_u(:) = 1
    print *,"Selecting reactions..."
    do i=1,nrea
       !reactants+products indexes
       rp = (/arr_r1(i), arr_r2(i), arr_r3(i), &
            arr_p1(i), arr_p2(i), arr_p3(i), arr_p3(i)/)
       !loop on reactants+products
       do j=1,size(rp)
          if(rp(j)==idx_dummy) cycle
          if(uspec(rp(j))==0) then
             arr_u(i) = 0
             exit
          end if
       end do
    end do
    print *,"Reaction used:",sum(arr_u)
    print *,"Reaction totl:",nrea

  end subroutine krome_authub
#ENDIFKROME

  !**************************
  function fex_check(n,Tgas)
    use krome_commons
    use krome_tabs
    implicit none
    integer::i
#KROME_rvars
    real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe_tab(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
#KROME_arrs
#KROME_arr_flux
       rrmax = max(rrmax, arr_flux(i))
    end do
    fex_check = rrmax
    
  end function fex_check

#IFKROME_useFlux
  !*************************
  subroutine flux_reduction(rrmax,threshold)
    use krome_commons
    implicit none
    integer::i
    real*8::rrmax,thold
    real*8,optional::threshold

    thold = 1d-8
    if(present(threshold)) thold = threshold
    
    do i=1,nrea
       if(arr_flux(i)>rrmax*thold) then
          arr_u(i) = 1
       else
          arr_u(i) = 0
       end if
    end do
 
  end subroutine flux_reduction
#ENDIFKROME

end module krome_reduction
