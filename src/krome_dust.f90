module krome_dust
contains
#IFKROME_useDust
  subroutine krome_init_dust(xdust,adust,ntot,alow_arg,aup_arg,phi_arg)
    !krome_init_dust: initialize the dust ditribution (xdust)
    ! and the dust bin mean sizes (adust). Arguments are
    ! ntot(number_of_dust_types) the total abundance per 
    ! each bin size, alow_arg the size of the smallest
    ! bin size, aup_arg the largest, phi_arg the exponent
    ! of the MNR power law.
    use krome_commons
    use krome_subs
    implicit none
    real*8,optional::alow_arg,aup_arg,phi_arg
    real*8::iphi1,c,phi1,abin(ndust+1),mass(nspec),xdust(ndust)
    real*8::alow,aup,phi,ntot(ndustTypes),adust(ndust),nd
    integer::i,j,ilow,iup

#KROME_dustParnerIndex

    !default values
    alow = 5d-8 !lower size (cm)
    aup = 2.5d-5 !upper size (cm)
    phi = -3.5d0 !MNR distribution exponent (with its sign)
    if(present(alow_arg)) alow = alow_arg
    if(present(aup_arg)) aup = alow_arg
    if(present(phi_arg)) phi = phi_arg
    mass(:) = get_mass()
    nd = ndust/ndustTypes
    !loop on dust types
    do j=1,ndustTypes
       ilow = nd * (j - 1) + 1 !lower index
       iup = nd * j !upper index
       phi1 = phi + 1.d0 !phi+1
       iphi1 = 1.d0/phi1 !1/phi
       
       !estimate normalization constant
       c = (phi1)/(aup**(phi1)-alow**(phi1))
       abin(1) = alow !set lower limit for the first bin
       !evaluates bin limits
       do i=2,nd+1
          abin(i) = (phi1/nd/c + abin(i-1)**phi1)**iphi1
       end do
       !evaluate means size of dust in each bin
       do i=1,nd
          adust(i+ilow-1) = abs(abin(i)-abin(i+1))*0.5d0
       end do
       !compute and normalize bin distribution
       xdust(ilow:iup) = c*adust(ilow:iup)**phi
       xdust(ilow:iup) = xdust(ilow:iup)/sum(xdust(ilow:iup)) * ntot(j)
       
       !evaluate dust-parnter ratio (e.g. 1dust=1e2 C atoms)
       krome_dust_partner_ratio(ilow:iup) = adust(ilow:iup)**3 &
            * 2.3d0  / mass(krome_dust_partner_idx(j))
       krome_dust_partner_ratio_inv(ilow:iup) = 1.d0 &
            / krome_dust_partner_ratio(ilow:iup)
    end do

    !compute the mass of the dust partner
    do j=1,ndustTypes
       krome_dust_partner_mass(j) =  mass(krome_dust_partner_idx(j))
    end do
    print *,"Dust initialized!"

  end subroutine krome_init_dust
#ENDIFKROME
end module krome_dust
