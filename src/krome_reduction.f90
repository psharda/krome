module krome_reduction
contains

#IFKROME_useTopology

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
    !this routine is not supported, use at your own risk
    ! the method employed here is different from Grassi+2012
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
