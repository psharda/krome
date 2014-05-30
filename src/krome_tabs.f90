module krome_tabs
contains

#IFKROME_useTabs
  subroutine make_ktab()
    !build the tabs from coefficients
    use krome_commons
    use krome_subs
    use krome_photo
    integer::i,j,ierror,kwarnup(nrea),kwarndown(nrea),pblock
    real*8::kk(nrea),valmax,n(nspec)
    logical::is_rank_zero
    !temperature limits
#KROME_logTlow
#KROME_logTup

    is_rank_zero = (krome_mpi_rank==0)

    !loop to create tabs (it may take a while)
    valmax = 1d0
    ierror = 0 !error count
    pblock = ktab_n/10 !ouput frequency
    if(is_rank_zero) print *,"KROME: creating tabs..."
    kwarnup(:) = 0 !store warnings
    kwarndown(:) = 0 !store warnings
    !loop on temperatures
    do i=1,ktab_n
       if(mod(i,pblock)==0 .and. is_rank_zero) print *,i/pblock*10,"%"
       ktab_T(i) = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)&
            +ktab_logTlow)
       n(:) = 0.d0
       n(idx_Tgas) = ktab_T(i)
       kk(:) = coe(n(:))
       !check for errors or discrepancies
       if((maxval(kk)>valmax.or.minval(kk)<0.d0)) then
          ierror = ierror + 1
          if(ierror==1.and. is_rank_zero) print '(a16,a5,2a11)',"",&
               "idx","Tgas","rate"
          do j=1,nrea
             if(kk(j)>valmax .and. kwarnup(j)==0) then
                kwarnup(j) = 1
                if(is_rank_zero) print '(a16,I5,2E11.3)', "WARNING: k>1.",&
                     j,ktab_T(i),kk(j)
             end if
             if(kk(j)<0.d0 .and. kwarndown(j)==0) then
                kwarndown(j) = 1
                if(is_rank_zero) print *,"WARNING: k<0.d0",j,ktab_T(i),kk(j)
             end if
          end do
       end if
       ktab(:,i) = kk(:)
    end do

    !store inverse values of deltaT to speed-up at runtime
    do i=1,ktab_n-1
       inv_ktab_T(i) = 1.d0 / (ktab_T(i+1)-ktab_T(i))
    end do

    !store inverse to go fast at runtime
    inv_ktab_idx = 1.d0 / (ktab_logTup - ktab_logTlow) * (ktab_n - 1)

  end subroutine make_ktab

  !************************
  subroutine check_tabs()
    use krome_commons
    use krome_subs
    integer::i,j,pblock,ii
    real*8::kk(nrea),kktab(nrea),Tgas,kmax,n(nspec),kold(nrea),dk
    logical::is_rank_zero

    is_rank_zero = (krome_mpi_rank==0)
    
    pblock = ktab_n/10 !write % every 10
    if(is_rank_zero) print *,"KROME: checking tabs..."
    !loop on tabs
    do i=1,ktab_n
       if(mod(i,pblock)==0.and.is_rank_zero) print *,i/pblock*10,"%" !output
       Tgas = 1d1**((i-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)+ktab_logTlow) 
       n(:) = 0.d0 !rates do not depends on densities
       n(idx_Tgas) = Tgas !rates depend on temperature
       kk(:) = coe(n(:)) !get rates
       kktab(:) = coe_tab(n(:)) !get rates from tabs
       kold(:) = 0.d0 !old rate value to skip discontinuities
       !loop on reactions
       do j=1,nrea
          kmax = kk(j)
          if(kmax>0.d0.and.kk(j)>0.d0) then
             dk = abs(kk(j)-kold(j))/(kold(j)+1d-40)
             if(abs(kk(j)-kktab(j))/kmax>1d-1.and.kmax>1d-12.and.dk<1d-1) then
                if(is_rank_zero) then
                   print *,"ERROR: wrong rate tables"
                   print *,"Rate index:",j
                   print *,"Temperature:",Tgas
                   print *,"Rate values:",kk(j),kktab(j)
                   print *,"Error:",abs(kk(j)-kktab(j))/kmax,&
                        "(must be close to zero)"

                   !dump graph
                   open(93,file="KROME_TAB_DUMP.dat",status="replace")
                   do ii=1,ktab_n
                      Tgas = 1d1**((ii-1)*(ktab_logTup-ktab_logTlow)/(ktab_n-1)&
                           +ktab_logTlow)
                      n(idx_Tgas) = Tgas !rates depend on temperature
                      kk(:) = coe(n(:))
                      kktab(:) = coe_tab(n(:))
                      write(93,'(99E12.3e3)') Tgas,kk(j),kktab(j)
                   end do
                   close(93)
                   print *,"Graph dump to KROME_TAB_DUMP.dat"
                   print *,"gnuplot command:"
                   print *," plot 'KROME_TAB_DUMP.dat' w l, '' u 1:3"
                   stop
                end if
             end if
          end if
       end do
       kold(:) = kk(:)
    end do
    if(is_rank_zero) print *,"KROME: tabs are ok!"

  end subroutine check_tabs
#ENDIFKROME

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_commons
    use krome_user_commons
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small,t
#KROME_define_vars

    Tgas = max(2.73d0,n(idx_Tgas))
    t = Tgas
    small = 0d0

#KROME_Tshortcuts
#KROME_init_vars

#IFKROME_useCustomCoe
    coe_tab(:) = #KROMEREPLACE_customCoeFunction
#ENDIFKROME

#IFKROME_useStandardCoe
    coe_tab(:) = coe(n(:))
#ENDIFKROME

#IFKROME_useTabs
    idx = (log10(Tgas)-ktab_logTlow) * inv_ktab_idx + 1
    idx = max(idx,1)
    idx = min(idx,ktab_n-1)
    coe_tab(:) = 0.d0
    do j=1,nrea
       coe_tab(j) = (Tgas-ktab_T(idx)) * inv_ktab_T(idx) *&
            (ktab(j,idx+1)-ktab(j,idx)) + ktab(j,idx)
    end do

    #KROME_noTabReactions
#ENDIFKROME

  end function coe_tab

end module krome_tabs
