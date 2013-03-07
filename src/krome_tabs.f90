module krome_tabs
contains
  subroutine make_ktab()
    use krome_commons
    use krome_subs
    integer::i,j,ierror,kwarnup(nrea),kwarndown(nrea),pblock
    real*8::kk(nrea),valmax,n(nspec)
    
    #KROME_logTlow
    #KROME_logTup

    valmax = 1d0
    ierror = 0
    pblock = ktab_n/10
    print *,"KROME: creating tabs..."
    kwarnup(:) = 0
    kwarndown(:) = 0
    do i=1,ktab_n
       if(mod(i,pblock)==0) print *,i/pblock*10,"%"
       ktab_T(i) = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)&
            +ktab_logTlow)
       n(:) = 0.d0
       n(idx_Tgas) = ktab_T(i)
       kk(:) = coe(n(:))
       if((maxval(kk)>valmax.or.minval(kk)<0.d0)) then
          ierror = ierror + 1
          if(ierror==1) print '(a16,a5,2a11)',"","idx","Tgas","rate"
          do j=1,nrea
             if(kk(j)>valmax .and. kwarnup(j)==0) then
                kwarnup(j) = 1
                print '(a16,I5,2E11.3)', "WARNING: k>1.",j,ktab_T(i),kk(j)
             end if
             if(kk(j)<0.d0 .and. kwarndown(j)==0) then
                kwarndown(j) = 1
                print *,"WARNING: k<0.d0",j,ktab_T(i),kk(j)
             end if
          end do
       end if
       ktab(:,i) = kk(:)
    end do

    do i=1,ktab_n-1
       inv_ktab_T(i) = 1.d0 / (ktab_T(i+1)-ktab_T(i))
    end do

    inv_ktab_idx = 1.d0 / (ktab_logTup - ktab_logTlow) * ktab_n
    
  end subroutine make_ktab
  
  !************************
  subroutine check_tabs()
    use krome_commons
    use krome_subs
    integer::i,j,pblock
    real*8::kk(nrea),kktab(nrea),Tgas,kmax,n(nspec)
    pblock = ktab_n/10
    print *,"KROME: checking tabs..."
    do i=1,ktab_n
       if(mod(i,pblock) == 0) print *,i/pblock*10,"%"
       Tgas = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)+ktab_logTlow) 
       n(:) = 0.d0
       n(idx_Tgas) = Tgas
       kk(:) = coe(n(:))
       kktab(:) = coe_tab(n(:))
       do j=1,nrea
          kmax = kk(j)
          if(kmax>0.d0.and.kk(j)>0.d0) then
             if(abs(kk(j)-kktab(j))/kmax>1d-2.and.kmax>1d-12) then
                print *,"ERROR: wrong rate tables"
                print *,"Rate index:",j
                print *,"Temperature:",Tgas
                print *,"Rate values:",kk(j),kktab(j)
                print *,"Error:",abs(kk(j)-kktab(j))/kmax,&
                     "(must be close to zero)"
                stop
             end if
          end if
       end do
    end do
    print *,"KROME: tabs are ok!"
    
  end subroutine check_tabs

  !***********************+
  function coe_tab(n)
    use krome_subs
    use krome_commons
    use krome_user_commons
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec)

    Tgas = max(2.73d0,n(idx_Tgas))

#IFKROME_useCustomCoe
    coe_tab(:) = #KROMEREPLACE_customCoeFunction
#ENDIFKROME

#IFKROME_useStandardCoe
    coe_tab(:) = coe(n(:))
#ENDIFKROME
    
#IFKROME_useTabs
    idx = (log10(Tgas)-ktab_logTlow) * inv_ktab_idx
    idx = max(idx,1)
    idx = min(idx,ktab_n-1)
    coe_tab(:) = 0.d0
    do j=1,nrea
       coe_tab(j) = (Tgas-ktab_T(idx)) * inv_ktab_T(idx) *&
            (ktab(j,idx+1)-ktab(j,idx)) + ktab(j,idx+1)
    end do
#ENDIFKROME

  end function coe_tab

end module krome_tabs
