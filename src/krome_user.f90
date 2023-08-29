module krome_user
#IFKROME_useBindC
  use iso_c_binding
#ENDIFKROME_useBindC
  implicit none

#KROME_header

#KROME_species

#KROME_cool_index

#KROME_heat_index

#KROME_common_alias

#KROME_constant_list

contains

#KROME_user_commons_functions

#IFKROME_customFex
  !*******************
  !custom FEX prototype function
  !subroutine myCustomFex(neq,tt,nin,dn)
  !  implicit none
  !  integer::neq
  !  real*8::tt,dn(neq),nin(neq)
  !
  !end subroutine myCustomFex
#ENDIFKROME

  !************************
  !returns the Tdust averaged over the number density
  ! as computed in the tables
  function krome_get_table_Tdust(x,Tgas) #KROME_bindC
    use krome_commons
    use krome_grfuncs
    implicit none
    #KROME_double_value :: Tgas
    #KROME_double :: x(nmols), krome_get_table_Tdust
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_table_Tdust = get_table_Tdust(n(:))

  end function krome_get_table_Tdust

#IFKROME_useBindC
  !MOCASSIN interface not used when C interface is active

#ELSEKROME_useBindC
  !**********************
  !convert from MOCASSIN abundances to KROME
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  ! imap: matrix position index map, integer
  ! returns KROME abundances (cm-3, real*8)
  function krome_convert_xmoc(xmoc,imap) result(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*4,intent(in):: xmoc(:,:)
    real*8::x(nmols),n(nspec)
    integer,intent(in)::imap(:)

    x(:) = 0d0

#KROME_xmoc_map

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
#IFKROME_has_electrons
    x(idx_e) = get_electrons(n(:))
#ENDIFKROME

  end function krome_convert_xmoc

  !*************************
  !convert from KROME abundances to MOCASSIN
  ! x: KROME abuances (cm-3, real*8)
  ! imap: matrix position index map, integer
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  subroutine krome_return_xmoc(x,imap,xmoc)
    use krome_commons
    implicit none
    real*8,intent(in)::x(nmols)
    real*4,intent(out)::xmoc(:,:)
    integer,intent(in)::imap(:)

    xmoc(:,:) = 0d0

#KROME_xmoc_map_return

  end subroutine krome_return_xmoc
#ENDIFKROME_useBindC

  !**********************
  !convert number density (cm-3) into column
  ! density (cm-2) using the specific density
  ! column method (see help for option
  ! -columnDensityMethod)
  ! num is the number density, x(:) is the species
  ! array, Tgas is the gas temperature
  ! If the method is not JEANS, x(:) and Tgas
  ! are dummy variables
  function krome_num2col(num,x,Tgas) #KROME_bindC
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    #KROME_double :: x(nmols),krome_num2col
    #KROME_double_value :: Tgas,num
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_num2col = num2col(num,n(:))

  end function krome_num2col

  !***********************
  !print on screen the current values of all phys variables
  subroutine krome_print_phys_variables() #KROME_bindC
    use krome_commons
    implicit none

#KROME_print_phys_variables

  end subroutine krome_print_phys_variables

#KROME_set_get_phys_functions


#IFKROME_useXrays
  !****************************
  !set value of J21xrays for tabulated
  ! heating and corresponding rate
  subroutine krome_set_J21xray(xarg) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double_value :: xarg

    J21xray = xarg

  end subroutine krome_set_J21xray
#ENDIFKROME

#KROME_cooling_functions

#IFKROME_use_coolingZ
  !***************************
  !dump the population of the Z cooling levels
  ! in the nfile file unit, using xvar as
  ! independent variable. alias of
  ! dump_cooling_pop subroutine
  subroutine krome_popcool_dump(xvar,nfile) #KROME_bindC
    use krome_cooling
    implicit none
    #KROME_double_value :: xvar
    #KROME_integer_value :: nfile

    call dump_cooling_pop(xvar,nfile)

  end subroutine krome_popcool_dump

#ENDIFKROME

#IFKROME_useTabsTdust

  !************************************
  ! when using 3D variables set the value of the 3rd dimension
  ! (usually Av). The other 2 dimensions are ntot and Tgas
  ! (take from internal KROME)
  subroutine krome_set_dust_table_3D_variable(variable_value)
    use krome_commons
    implicit none
    real*8,intent(in)::variable_value

    dust_table_AvVariable_log = log10(variable_value)

  end subroutine krome_set_dust_table_3D_variable

  !*****************************
  !get averaged Tdust from tables, with x(:) species array
  ! of size krome_nmols, and Tgas the gas temperature
  function krome_get_Tdust(x,Tgas) #KROME_bindC
    use krome_commons
    use krome_fit
    implicit none
    #KROME_double :: x(nmols), krome_get_Tdust
    #KROME_double_value :: Tgas
    real*8 :: ntot

    ntot = sum(x(1:nmols))
#IFKROME_dust_table_2D
    krome_get_Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), dust_tab_Tgas(:), &
         dust_tab_Tdust(:,:), dust_mult_ngas, dust_mult_Tgas, &
         log10(ntot), log10(Tgas))
#ENDIFKROME_dust_table_2D

#IFKROME_dust_table_3D
    krome_get_Tdust = 1d1**fit_anytab3D(dust_tab_ngas(:), dust_tab_Tgas(:), &
         dust_tab_AvVariable(:), dust_tab_Tdust(:,:,:), dust_mult_ngas, &
         dust_mult_Tgas, dust_mult_AvVariable, &
         log10(ntot), log10(Tgas), dust_table_AvVariable_log)
#ENDIFKROME_dust_table_3D

  end function krome_get_Tdust

#ENDIFKROME

#IFKROME_useDust

  !*************************
  !this subroutine sets the dust distribution in the range
  ! alow_arg, aup_arg, using power law with exponent phi_arg.
  ! All these arguments are optional, execept for x(:) of size
  ! krome_nmols that represents the number densitites of the
  ! chemical species, and dust_gas_ratio. Check the code below
  ! for defaults
  subroutine krome_init_dust_distribution(x,dust_gas_ratio,alow_arg,&
       aup_arg,phi_arg) #KROME_bindC
    use krome_dust
    use krome_commons
    #KROME_double_value_optional :: alow_arg,aup_arg,phi_arg
    #KROME_double_value :: dust_gas_ratio
    #KROME_double :: x(nmols)
    real*8::alow,aup,phi

    !default values
#IFKROME_useBindC
    alow = alow_arg
    aup  = aup_arg
    phi  = phi_arg
#ELSEKROME_useBindC
    alow = 5d-7 !lower size (cm)
    aup = 2.5d-5 !upper size (cm)
    phi = -3.5d0 !MNR distribution exponent (with its sign)

    if(present(alow_arg)) alow = alow_arg
    if(present(aup_arg)) aup = aup_arg
    if(present(phi_arg)) phi = phi_arg
#ENDIFKROME_useBindC

    call set_dust_distribution(x(:),dust_gas_ratio,alow,aup,phi)

  end subroutine krome_init_dust_distribution

  !*****************************
  !this function returns an array of size krome_ndust
  ! that contains the amount of dust per bin in 1/cm3.
  function krome_get_dust_distribution()
    use krome_commons
    implicit none
    real*8 :: krome_get_dust_distribution(ndust)

    krome_get_dust_distribution(:) = xdust(:)

  end function krome_get_dust_distribution

  !*****************************
  !this function sets the dust distribution with an array
  ! that contains the amount of dust per bin in 1/cm3.
  ! each array element represents a dust bin.
  subroutine krome_set_dust_distribution(arg) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double :: arg(ndust)

    xdust(:) = arg(:)

  end subroutine krome_set_dust_distribution

  !******************************
  !this function returns an array of size krome_ndust
  ! that contains the size of the dust bins in cm
  function krome_get_dust_size()
    use krome_commons
    implicit none
    real*8 :: krome_get_dust_size(ndust)

    krome_get_dust_size(:) = krome_dust_asize(:)

  end function krome_get_dust_size

  !******************************
  !this function set the sizes of the dust with an array
  ! of size krome_ndust that contains the size of the
  ! dust bins in cm
  subroutine krome_set_dust_size(arg) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double :: arg(ndust)

    krome_dust_asize(:) = arg(:)
    krome_dust_asize2(:) = arg(:)**2
    krome_dust_asize3(:) = arg(:)**3

  end subroutine krome_set_dust_size

  !************************
  !this function sets the default temperature (arg)
  ! for all the dust bins, K
  subroutine krome_set_Tdust(arg) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double_value :: arg

    krome_dust_T(:) = arg

  end subroutine krome_set_Tdust

  !************************
  !this function sets the temperature
  ! for all the dust bins but using an array
  ! of size krome_ndust, K
  subroutine krome_set_Tdust_array(arr) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double :: arr(ndust)

    krome_dust_T(:) = arr(:)

  end subroutine krome_set_Tdust_array

  !*********************
  !returns the Tdust averaged over the dust number density
  function krome_get_averaged_Tdust() #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double :: krome_get_averaged_Tdust

    krome_get_averaged_Tdust = sum(xdust(:)*krome_dust_T(:))/sum(xdust(:))

  end function krome_get_averaged_Tdust

  !****************************
  ! scales the dust distribution by multiplying it by the
  ! real*8 value xscale
  subroutine krome_scale_dust_distribution(xscale) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_double_value :: xscale

    xdust(:) = xdust(:) * xscale

  end subroutine krome_scale_dust_distribution

  !***********************
  !returns an array of size krome_ndust containing the
  ! dust temperatures, K. Each array element is a dust bin
  function krome_get_Tdust()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tdust(ndust)

    krome_get_Tdust(:) = krome_dust_T(:)

  end function krome_get_Tdust

  !***************************
  !this subroutine sets as a scalar (xarg) all the
  ! surface species for the given idx_base in the
  ! species array x:
  ! e.g. krome_set_surface(x(:),1d3,krome_idx_OH_dust)
  subroutine krome_set_surface(x,xarg,idx_base) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_integer_value :: idx_base
    #KROME_double_value :: xarg
    #KROME_double :: x(nmols)

    x(idx_base:idx_base+ndust-1) = xarg

  end subroutine krome_set_surface

  !***************************
  !this subroutine sets as a scalar (xarg) all the
  ! surface species for the given idx_base in the
  ! species array x, normalized by the amount of dust
  ! in each bin:
  ! e.g. krome_set_surface_norm(x(:),1d3,krome_idx_OH_dust)
  subroutine krome_set_surface_norm(x,xarg,idx_base) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_integer_value :: idx_base
    #KROME_double_value :: xarg
    #KROME_double :: x(nmols)

    x(idx_base:idx_base+ndust-1) = xarg*xdust(:)/sum(xdust(:))

  end subroutine krome_set_surface_norm

  !***************************
  !this subroutine sets as a vector (xarr) all the
  ! surface species for the given idx_base in the
  ! species array x. The size of the array xarr is ndust.
  ! e.g. krome_set_surface_array(x(:),arr(:),krome_idx_OH_dust)
  subroutine krome_set_surface_array(x,xarr,idx_base) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_integer_value :: idx_base
    #KROME_double :: xarr(ndust),x(nmols)

    x(idx_base:idx_base+ndust-1) = xarr(:)

  end subroutine krome_set_surface_array

  !***************************
  !this function gets the total amount of surface
  ! species for the given idx_base in the
  ! species array x.
  ! e.g. xx = krome_get_surface(x(:),krome_idx_OH_dust)
  function krome_get_surface(x,idx_base) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_integer_value :: idx_base
    #KROME_double :: x(nmols), krome_get_surface

    krome_get_surface = sum(x(idx_base:idx_base+ndust-1))

  end function krome_get_surface
#ENDIFKROME

  !*****************************
  !dump the data for restart (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_store(x,Tgas,dt) #KROME_bindC
    use krome_commons
    implicit none
    integer::nfile,i
    #KROME_double :: x(nmols)
    #KROME_double_value :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="replace")
    !dump temperature
    write(nfile,*) Tgas
    write(nfile,*) dt
    !dump species
    do i=1,nmols
       write(nfile,*) x(i)
    end do
#IFKROME_useDust
    !dump dust
    do i=1,ndust
       write(nfile,*) xdust(i)
    end do
#ENDIFKROME
    close(nfile)

  end subroutine krome_store

  !*****************************
  !restore the data from a dump (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_restore(x,Tgas,dt) #KROME_bindC
    use krome_commons
    implicit none
    integer::nfile,i
    #KROME_double :: x(nmols)
    #KROME_double_value :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="old")
    !restore temperature
    read(nfile,*) Tgas
    read(nfile,*) dt
    !restore species
    do i=1,nmols
       read(nfile,*) x(i)
    end do
#IFKROME_useDust
    !restore dust
    do i=1,ndust
       read(nfile,*) xdust(i)
    end do
#ENDIFKROME
    close(nfile)

  end subroutine krome_restore

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on() #KROME_bindC
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off() #KROME_bindC
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

#IFKROME_usePhotoBins
  !************************
  ! prepares tables for cross sections and
  ! photorates
  subroutine krome_calc_photobins() #KROME_bindC
    use krome_photo
    call calc_photobins()
  end subroutine krome_calc_photobins

  !****************************
  ! set the energy per photo bin
  ! eV/cm2/sr
  subroutine krome_set_photoBinJ(phbin) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
    #KROME_double :: phbin(nPhotoBins)
    photoBinJ(:) = phbin(:)
    photoBinJ_org(:) = phbin(:) !for restore

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBinJ

  !*************************
  ! set the energy (frequency) of the photobin
  ! as left-right limits in eV
  subroutine krome_set_photobinE_lr(phbinleft,phbinright,Tgas) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
#KROME_double :: phbinleft(nPhotoBins),phbinright(nPhotoBins)
#KROME_double_value_optional::Tgas
    real*8::bTgas

#IFKROME_useBindC
    bTgas = tgas
#ELSEKROME_useBindC
    !default Tgas for broadening
    bTgas = 1d1
    if(present(Tgas)) then
       bTgas = Tgas
    end if
#ENDIFKROME_useBindC

    !$omp parallel
    photoBinEleft(:) = phbinleft(:)
    photoBinEright(:) = phbinright(:)
    photoBinEmid(:) = 0.5d0*(phbinleft(:)+phbinright(:))
    photoBinEdelta(:) = phbinright(:)-phbinleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)
    !$omp end parallel

    !initialize xsecs table
    call init_photoBins(bTgas)

  end subroutine krome_set_photobinE_lr

  !*************************
  ! set the energy (frequency) of photobins
  ! when contiguous. Left and right limits are automatically
  ! extracted. Energy in eV
  subroutine krome_set_photobinE_limits(phbinLimits,Tgas) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
    #KROME_double :: phbinLimits(nPhotoBins+1)
    #KROME_double_value_optional::Tgas
    real*8::phl(nPhotoBins),phr(nPhotoBins),bTgas

#IFKROME_useBindC
    bTgas = tgas
#ELSEKROME_useBindC
    !default Tgas for broadening
    bTgas = 1d1
    if(present(Tgas)) then
       bTgas = Tgas
    end if
#ENDIFKROME_useBindC
    phl(:) = phbinLimits(1:nPhotoBins)
    phr(:) = phbinLimits(2:nPhotoBins+1)

    call krome_set_photobinE_lr(phl(:),phr(:),bTgas)

  end subroutine krome_set_photobinE_limits

  !*******************************
  !set the energy (eV) of the photobin according
  ! to MOCASSIN way (position and width array)
  subroutine krome_set_photobinE_moc(binPos,binWidth,Tgas) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
#KROME_double :: binPos(nPhotoBins),binWidth(nPhotoBins)
#KROME_double_value_optional::Tgas
    real*8::bTgas

#IFKROME_useBindC
    bTgas = Tgas
#ELSEKROME_useBindC
    bTgas = 1d1
    if(present(Tgas)) then
       bTgas = Tgas
    end if
#ENDIFKROME_useBindC

    !$omp parallel
    photoBinEleft(:) = binPos(:)-binWidth(:)/2d0
    photoBinEright(:) = binPos(:)+binWidth(:)/2d0
    photoBinEmid(:) = binPos(:)
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)
    !$omp end parallel

    !initialize xsecs table
    call init_photoBins(bTgas)

  end subroutine krome_set_photobinE_moc

  !********************************
  ! set the energy (eV) of the photobin
  ! linearly from lowest to highest energy value
  ! in eV
  subroutine krome_set_photobinE_lin(lower,upper,Tgas) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
#KROME_double_value :: lower,upper
#KROME_double_value_optional::Tgas
    real*8::dE,bTgas
    integer::i

#IFKROME_useBindC
    bTgas = Tgas
#ELSEKROME_useBindC
    bTgas = 1d1
    if(present(Tgas)) then
       bTgas = Tgas
    end if
#ENDIFKROME_useBindC

    !$omp parallel
    dE = abs(upper-lower)/nPhotoBins
    !$omp end parallel
    do i=1,nPhotoBins
       !$omp parallel
       photoBinEleft(i) = dE*(i-1) + lower
       photoBinEright(i) = dE*i + lower
       photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
       !$omp end parallel
    end do
    !$omp parallel
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)
    !$omp end parallel

    !initialize xsecs table
    call init_photoBins(bTgas)

  end subroutine krome_set_photobinE_lin

  !********************************
  ! set the energy (eV) of the photobin
  ! logarithmically from lowest to highest energy value
  ! in eV
  subroutine krome_set_photobinE_log(lower,upper,Tgas) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
#KROME_double_value :: lower,upper
#KROME_double_value_optional::Tgas
    real*8::dE,logup,loglow,bTgas
    integer::i

#IFKROME_useBindC
    bTgas = Tgas
#ELSEKROME_useBindC
    bTgas = 1d1
    if(present(Tgas)) then
       bTgas = Tgas
    end if
#ENDIFKROME_useBindC

    if(lower.ge.upper) then
       print *,"ERROR: in  krome_set_photobinE_log lower >= upper limit!"
       stop
    end if
    loglow = log10(lower)
    logup = log10(upper)
    !$omp parallel
    dE = 1d1**(abs(logup-loglow)/nPhotoBins)
    !$omp end parallel
    do i=1,nPhotoBins
       !$omp parallel
       photoBinEleft(i) = 1d1**((i-1)*(logup-loglow)/nPhotoBins + loglow)
       photoBinEright(i) = 1d1**(i*(logup-loglow)/nPhotoBins + loglow)
       photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
       !$omp end parallel
    end do
    !$omp parallel
    photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
    photoBinEidelta(:) = 1d0/photoBinEdelta(:)
    !$omp end parallel

    !initialize xsecs table
    call init_photoBins(bTgas)

  end subroutine krome_set_photobinE_log

  !*********************************
  !returns an array containing the flux for each photo bin
  ! in eV/cm2/sr
  function krome_get_photoBinJ()
    use krome_commons
    real*8 :: krome_get_photoBinJ(nPhotoBins)
    krome_get_photoBinJ(:) = photoBinJ(:)
  end function krome_get_photoBinJ

  !*********************************
  !get an array containing all the left positions
  ! of the photobins, eV
  function krome_get_photoBinE_left()
    !returns an array of size krome_nPhotoBins with the
    ! left energy limits (eV)
    use krome_commons
    real*8 :: krome_get_photoBinE_left(nPhotoBins)
    krome_get_photoBinE_left(:) = photoBinEleft(:)
  end function krome_get_photoBinE_left

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! right energy limits (eV)
  function krome_get_photoBinE_right()
    use krome_commons
    real*8 :: krome_get_photoBinE_right(nPhotoBins)
    krome_get_photoBinE_right(:) = photoBinEright(:)
  end function krome_get_photoBinE_right

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! middle energy values (eV)
  function krome_get_photoBinE_mid()
    use krome_commons
    real*8 :: krome_get_photoBinE_mid(nPhotoBins)
    krome_get_photoBinE_mid(:) = photoBinEmid(:)
  end function krome_get_photoBinE_mid

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! bin span (eV)
  function krome_get_photoBinE_delta()
    use krome_commons
    real*8 :: krome_get_photoBinE_delta(nPhotoBins)
    krome_get_photoBinE_delta(:) = photoBinEdelta(:)
  end function krome_get_photoBinE_delta

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! inverse of the bin span (1/eV)
  function krome_get_photoBinE_idelta()
    use krome_commons
    real*8 :: krome_get_photoBinE_idelta(nPhotoBins)
    krome_get_photoBinE_idelta(:) = photoBinEidelta(:)
  end function krome_get_photoBinE_idelta

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! integrated photo rates (1/s)
  function krome_get_photoBin_rates()
    use krome_commons
    real*8 :: krome_get_photoBin_rates(nPhotoRea)
    krome_get_photoBin_rates(:) = photoBinRates(:)
  end function krome_get_photoBin_rates

  !*********************************
  !returns an array of size krome_nPhotoBins containing
  ! the cross section (cm2) of the idx-th photoreaction
  function krome_get_xsec(idx)
    use krome_commons
    implicit none
    real*8 :: krome_get_xsec(nPhotoBins)
    integer :: idx

    krome_get_xsec(:) = photoBinJTab(idx,:)

  end function krome_get_xsec

  !*********************************
  !returns an array of size krome_nPhotoBins with the
  ! integrated photo heatings (erg/s)
  function krome_get_photoBin_heats()
    use krome_commons
    implicit none
    real*8 :: krome_get_photoBin_heats(nPhotoRea)
    krome_get_photoBin_heats(:) = photoBinHeats(:)

  end function krome_get_photoBin_heats

  !****************************
  !multiply all photobins by a factor real*8 xscale
  subroutine krome_photoBin_scale(xscale) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
    #KROME_double_value :: xscale

    photoBinJ(:) = photoBinJ(:) * xscale

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale

  !****************************
  !multiply all photobins by a real*8 array xscale(:)
  ! of size krome_nPhotoBins
  subroutine krome_photoBin_scale_array(xscale) #KROME_bindC
    use krome_commons
    use krome_photo
    implicit none
    #KROME_double :: xscale(nPhotoBins)

    photoBinJ(:) = photoBinJ(:) * xscale(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_photoBin_scale_array

  !********************************
  !restore the original flux (i.e. undo any rescale).
  ! the flux is automatically stored by the functions
  ! that set the flux, or by the function
  ! krome_photoBin_store()
  subroutine krome_photoBin_restore() #KROME_bindC
    use krome_commons
    implicit none

    photoBinJ(:) = photoBinJ_org(:)

  end subroutine krome_photoBin_restore

  !**********************
  !store flux to be restored with the subroutine
  ! krome_photoBin_restore later
  subroutine krome_photoBin_store() #KROME_bindC
    use krome_commons
    implicit none

    photoBinJ_org(:) = photoBinJ(:)

  end subroutine krome_photoBin_store

  !*********************
  !load flux radiation from a two-columns file
  ! energy/eV, flux/(eV/cm2/sr)
  ! Flux is interpolated over the existing binning
  ! constant-area method
  subroutine krome_load_photoBin_file_2col(fname, logarithmic)
    use krome_commons
    implicit none
    integer,parameter::imax=int(1e4)
    character(len=*) :: fname
    logical, optional :: logarithmic
    logical :: is_log
    integer::unit,ios,icount,j,i
    real*8::xtmp(imax),ftmp(imax),intA,eL,eR
    real*8::xL,xR,pL,pR,fL,fR,Jflux(nPhotoBins)
    real*8::a,b

    if(present(logarithmic)) then
      is_log = logarithmic
    else
      is_log = .false.
    end if

    !open file to read
    open(newunit=unit,file=trim(fname),iostat=ios)
    if(ios/=0) then
       print *,"ERROR: problems reading "//trim(fname)
       stop
    end if

    !read file line by line and store to temporary
    ftmp(:) = 0d0
    icount = 1
    do
       read(unit,*,iostat=ios) xtmp(icount), ftmp(icount)
       if(ios/=0) exit
       icount = icount + 1
    end do
    close(unit)
    icount = icount - 1

    if(is_log) ftmp = log(merge(ftmp,1d-40,ftmp>0d0))
    !loop on photobins for interpolation
    do j=1,nPhotoBins
       intA = 0d0
       !photobin limits
       eL = photoBinEleft(j)
       eR = photoBinEright(j)
       !loop on flux bins
       do i=1,icount-1
          !flux bin limits
          xL = xtmp(i)
          xR = xtmp(i+1)
          !if outside the bin skip
          if((xR<eL).or.(xL>eR)) cycle
          !get the interval limit (consider partial overlapping)
          pL = max(xL,eL)
          pR = min(xR,eR)
          if(is_log) then
            pL = log(max(pL,1d-40))
            pR = log(max(pR,1d-40))
            xL = log(max(xL,1d-40))
            xR = log(max(xR,1d-40))
          end if
          !interpolate to get the flux at the interval limit
          fL = (ftmp(i+1)-ftmp(i))*(pL-xL)/(xR-xL)+ftmp(i)
          fR = (ftmp(i+1)-ftmp(i))*(pR-xL)/(xR-xL)+ftmp(i)
          if(is_log) then
            fL = exp(fL)
            fR = exp(fR)
            pL = exp(pL)
            pR = exp(pR)
          end if

          !compute area of the overlapped area
          intA = intA + (fL+fR)*(pR-pL)/2d0
       end do
       !distribute the flux in the photobin
       Jflux(j) = intA / (eR-eL)
    end do

    !initialize intensity according to data
    call krome_set_photoBinJ(Jflux(:))

  end subroutine krome_load_photoBin_file_2col

  !********************************
  !load the radiation bins from the file fname
  ! data should be a 3-column file with
  ! energy Left (eV), energy Right (eV)
  ! intensity (eV/cm2/sr).
  ! This subroutine sets also the bin-size
  subroutine krome_load_photoBin_file(fname) !!#KROME_bindC !! not yet callable from C
    use krome_commons
    implicit none
    integer::ios,icount
    character(len=*) :: fname
    real*8::tmp_El(nPhotoBins),tmp_Er(nPhotoBins)
    real*8::rout(3),tmp_J(nPhotoBins)

    !open file and check for errors
    open(33,file=fname,status="old",iostat=ios)
    if(ios.ne.0) then
       print *,"ERROR: problem opening "//fname//"!"
       print *," (e.g. file not found)"
       stop
    end if

    icount = 0 !count valid line
    !loop on file
    do
       read(33,*,iostat=ios) rout(:)
       if(ios==-1) exit !EOF
       if(ios.ne.0) cycle !skip comments
       icount = icount + 1
       if(icount>nPhotoBins) exit !can't load more than nPhotoBins
       tmp_El(icount) = rout(1) !energy L eV
       tmp_Er(icount) = rout(2) !energy R eV
       !check if left interval is before right
       if(tmp_El(icount)>tmp_Er(icount)) then
          print *,"ERROR: in file "//fname//" left"
          print *, " interval larger than right one!"
          print *,tmp_El(icount),tmp_Er(icount)
          stop
       end if
       tmp_J(icount) = rout(3) !intensity eV/cm2/sr
    end do
    close(33)

    !file data lines should be the same number of the photobins
    if(icount/=nPhotoBins) then
       print *,"ERROR: the number of data lines in the file"
       print *," "//fname//" should be equal to the number of"
       print *," photobins ",nPhotoBins
       print *,"Found",icount
       stop
    end if

    !initialize interval and intensity according to data
#IFKROME_useBindC
    call krome_set_photobinE_lr(tmp_El(:),tmp_Er(:),1d1)
#ELSEKROME_useBindC
    call krome_set_photobinE_lr(tmp_El(:),tmp_Er(:))
#ENDIFKROME_useBindC
    call krome_set_photoBinJ(tmp_J(:))

  end subroutine krome_load_photoBin_file

  !**********************************
  !this subroutine sets an Hardt+Madau flux in the
  ! energy limits lower_in, upper_in (eV, log-spaced)
  subroutine krome_set_photoBin_HMlog(lower_in,upper_in) #KROME_bindC
    use krome_commons
    use krome_photo
    use krome_subs
    use krome_fit
    implicit none
    real*8::z(59),energy(500),HM(59,500)
    real*8::z_mul,energy_mul,x,lower,upper
    real*8,parameter::limit_lower = 0.1237d0
    real*8,parameter::limit_upper = 4.997d7
    real*8,parameter::limit_redshift = 15.660d0
    #KROME_double_value_optional::lower_in,upper_in
    integer::i

#IFKROME_useBindC
    lower = lower_in
    upper = upper_in
#ELSEKROME_useBindC
    lower = limit_lower
    upper = limit_upper
    if(present(lower_in)) lower = lower_in
    if(present(upper_in)) upper = upper_in
#ENDIFKROME_useBindC

    if(phys_zredshift>limit_redshift) then
       print *,"ERROR: redshift out of range in HM"
       print *,"redshift:",phys_zredshift
       print *,"limit:",limit_redshift
       stop
    end if

    if(lower<limit_lower .or. upper>limit_upper) then
       print *,"ERROR: upper or lower limit out of range in HM."
       print *,"lower limit (eV):",limit_lower
       print *,"upper limit (eV):",limit_upper
       stop
    end if

#IFKROME_useBindC
    call krome_set_photoBinE_log(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_log(lower,upper)
#ENDIFKROME_useBindC

    call init_anytab2D("krome_HMflux.dat", z(:), energy(:), &
         HM(:,:), z_mul, energy_mul)

    do i=1,nPhotoBins
       x = log10(photoBinEmid(i)) !log(eV)
       photoBinJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
            z_mul, energy_mul, phys_zredshift, x)
    end do

    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_HMlog

  !**********************************
  !this subroutine ADD an Hardt+Madau flux to the current radiation
  ! in the energy limits lower_in, upper_in (eV), It assumes
  ! the current binning
  subroutine krome_set_photoBin_HMCustom(lower_in,upper_in,additive) #KROME_bindC
    use krome_commons
    use krome_photo
    use krome_subs
    use krome_fit
    implicit none
    real*8::z(59),energy(500),HM(59,500)
    real*8::z_mul,energy_mul,x,lower,upper
    real*8::photoTmpJ(nPhotoBins)
    real*8,parameter::limit_lower = 0.1237d0
    real*8,parameter::limit_upper = 4.997d7
    real*8,parameter::limit_redshift = 15.660d0
    #KROME_bool_optional::additive
    logical::add
    #KROME_double_value_optional :: lower_in,upper_in
    integer::i

#IFKROME_useBindC
    lower = lower_in
    upper = upper_in
#ELSEKROME_useBindC
    lower = limit_lower
    upper = limit_upper
    if(present(lower_in)) lower = lower_in
    if(present(upper_in)) upper = upper_in
#ENDIFKROME_useBindC

    add = .false.
    if(present(additive)) add = additive

    if(phys_zredshift>limit_redshift) then
       print *,"ERROR: redshift out of range in HM"
       print *,"redshift:",phys_zredshift
       print *,"limit:",limit_redshift
       stop
    end if

    if(lower<limit_lower .or. upper>limit_upper) then
       print *,"ERROR: upper or lower limit out of range in HM."
       print *,"lower limit (eV):",limit_lower
       print *,"upper limit (eV):",limit_upper
       stop
    end if

    call init_anytab2D("krome_HMflux.dat", z(:), energy(:), &
         HM(:,:), z_mul, energy_mul)

    do i=1,nPhotoBins
       x = log10(photoBinEmid(i)) !log(eV)
       photoTmpJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
            z_mul, energy_mul, phys_zredshift, x)
    end do

    !add flux to already-present flux if optional argument
    if(add) then
       photoBinJ(:) = photoBinJ(:) + photoTmpJ(:)
    else
       photoBinJ(:) = photoTmpJ(:)
    end if

    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_HMCustom

  !**********************************
  !set the flux as a black body with temperature Tbb (K)
  ! in the range lower to upper (eV),  linear-spaced
  subroutine krome_set_photoBin_BBlin(lower,upper,Tbb) #KROME_bindC
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_phfuncs
    implicit none
    #KROME_double_value :: lower,upper,Tbb
    real*8::x
    integer::i

#IFKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper)
#ENDIFKROME_useBindC

    !eV/cm2/sr
    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       photoBinJ(i) = planckBB(x,Tbb)
    end do
    photoBinJ_org(:) = photoBinJ(:)

    call calc_photobins()

  end subroutine krome_set_photoBin_BBlin

  !**********************************
  !set the flux as a black body with temperature Tbb (K)
  ! in the range lower to upper (eV), log-spaced
  subroutine krome_set_photoBin_BBlog(lower,upper,Tbb) #KROME_bindC
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_phfuncs
    implicit none
    #KROME_double_value :: lower,upper,Tbb
    real*8::x,xmax,xexp,Jlim
    integer::i

    !limit for the black body intensity to check limits
    Jlim = 1d-3

#IFKROME_useBindC
    call krome_set_photoBinE_log(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_log(lower,upper)
#ENDIFKROME_useBindC

    !eV/cm2/sr
    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       photoBinJ(i) = planckBB(x,Tbb)
    end do
    photoBinJ_org(:) = photoBinJ(:)

    !uncomment this below for additional control
!!$    !find the maximum using Wien's displacement law
!!$    xmax = Tbb/2.8977721d-1 * clight * planck_eV !eV
!!$
!!$    if(xmax<lower) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is below the lowest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"lowest (eV)",lower
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(xmax>upper) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is above the highest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"highest (eV)",upper
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(1)>Jlim) then
!!$       print *,"WARNING: lower bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(1)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(nPhotoBins)>Jlim) then
!!$       print *,"WARNING: upper bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(nPhotoBins)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_BBlog

  !*************************************
  !set the BB spectrum and the limits using bisection
  subroutine krome_set_photoBin_BBlog_auto(Tbb) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_constants
    use krome_phfuncs
    implicit none
    #KROME_double_value :: Tbb
    real*8::xlow,xup,eps,xmax,J0,J1,x0,x1,xm,Jm
    eps = 1d-6

    !Rayleigh-Jeans approximation for the minimum energy
    xlow = planck_eV*clight*sqrt(.5d0/Tbb/boltzmann_eV*eps)

    !find energy of the Wien maximum (eV)
    xmax = Tbb / 2.8977721d-1 * clight * planck_eV

    !bisection to find the maximum
    x0 = xmax
    x1 = 2.9d2*Tbb*boltzmann_eV
    J0 = planckBB(x0,Tbb) - eps
    J1 = planckBB(x1,Tbb) - eps
    if(J0<0d0.or.J1>0d0) then
       print *,"ERROR: problems with auto planck bisection!"
       stop
    end if

    do
       xm = 0.5d0*(x0+x1)
       Jm = planckBB(xm,Tbb) - eps
       if(Jm>0d0) x0 = xm
       if(Jm<0d0) x1 = xm
       if(abs(Jm)<eps*1d-3) exit
    end do
    xup = xm

    !initialize BB radiation using the values found
    call krome_set_photoBin_BBlog(xlow,xup,Tbb)

  end subroutine krome_set_photoBin_BBlog_auto

  !*********************************
  !return the ratio between the current flux an Draine's
  function krome_get_ratioFluxDraine()
    use krome_subs
    use krome_phfuncs
    implicit none
    real*8::krome_get_ratioFluxDraine

    krome_get_ratioFluxDraine = get_ratioFluxDraine()

  end function krome_get_ratioFluxDraine

  !**********************************
  !set the flux as Draine's function
  ! in the range lower to upper (eV). the spacing is linear
  subroutine krome_set_photoBin_draineLin(lower,upper) #KROME_bindC
    use krome_commons
    use krome_photo
    use krome_constants
    #KROME_double_value :: upper,lower
    real*8::x
    integer::i

#IFKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper)
#ENDIFKROME_useBindC

    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       !eV/cm2/sr
       if(x<13.6d0.and.x>5d0) then
          photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
               * x *planck_eV
       else
          photoBinJ(i) = 0d0
       end if
    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLin

  !**************************
  !set the flux as Draine's function
  ! in the range lower to upper (eV). log-spaced
  subroutine krome_set_photoBin_draineLog(lower,upper) #KROME_bindC
    use krome_commons
    use krome_photo
    use krome_constants
    #KROME_double_value :: upper,lower
    real*8::x
    integer::i

#IFKROME_useBindC
    call krome_set_photoBinE_log(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_log(lower,upper)
#ENDIFKROME_useBindC

    do i=1,nPhotoBins
       x = photoBinEmid(i) !eV
       !eV/cm2/sr/s/Hz
       if(x<13.6d0.and.x>5d0) then
          photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
               * x *planck_eV
       else
          photoBinJ(i) = 0d0
       end if
    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineLog

  !**************************
  !set the flux as Draine's function with the current binning
  ! Note: you have to set the binning first
  subroutine krome_set_photoBin_draineCustom() #KROME_bindC
    use krome_commons
    use krome_photo
    use krome_constants
    real*8::xL,xR,f1,f2
    integer::i

    !return error if binning is not set
    if(maxval(photoBinEmid)==0d0) then
       print *,"ERROR: not initialized binning in draineCustom!"
       stop
    end if

    !loop on bins
    do i=1,nPhotoBins
       !eV/cm2/sr
       if(photoBinEright(i)<=13.6d0.and.photoBinEleft(i)>=5d0) then
          xL = photoBinEleft(i) !eV
          xR = photoBinEright(i) !eV
       elseif(photoBinEleft(i)<5d0.and.photoBinEright(i)>5d0) then
          xL = 5d0 !eV
          xR = photoBinEright(i) !eV
       elseif(photoBinEleft(i)<13d0.and.photoBinEright(i)>13.6d0) then
          xL = photoBinEleft(i) !eV
          xR = 13.6d0 !eV
       else
          xL = 0d0
          xR = 0d0
       end if
       f1 = (1.658d6*xL - 2.152d5*xL**2 + 6.919d3*xL**3) &
            * planck_eV
       f2 = (1.658d6*xR - 2.152d5*xR**2 + 6.919d3*xR**3) &
            * planck_eV
       photoBinJ(i) = (f1+f2)/2d0 * (xR-xL)/(photoBinEright(i)-photoBinEleft(i)) * photoBinEmid(i)

    end do

    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_draineCustom

  !**************************
  !set the flux as power-law (J21-style)
  ! in the range lower to upper (eV). linear-spaced
  subroutine krome_set_photoBin_J21lin(lower,upper) #KROME_bindC
    use krome_commons
    use krome_photo
    #KROME_double_value :: upper,lower

#IFKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_lin(lower,upper)
#ENDIFKROME_useBindC

    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/sr
    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21lin

  !**************************
  !set the flux as power-law (J21-style)
  ! in the range lower to upper (eV). the spacing is logarithmic
  subroutine krome_set_photoBin_J21log(lower,upper) #KROME_bindC
    use krome_commons
    use krome_photo
    #KROME_double_value :: upper,lower

#IFKROME_useBindC
    call krome_set_photoBinE_log(lower,upper,1d1)
#ELSEKROME_useBindC
    call krome_set_photoBinE_log(lower,upper)
#ENDIFKROME_useBindC

    photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/sr
    photoBinJ_org(:) = photoBinJ(:)

    !compute rates
    call calc_photobins()

  end subroutine krome_set_photoBin_J21log

  !*****************************
  !get the opacity tau corresponding to x(:)
  ! chemical composition. The column density
  ! is computed using the expression in the
  ! num2col(x) function.
  ! An array of size krome_nPhotoBins is returned
  function krome_get_opacity(x,Tgas)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols),krome_get_opacity(nPhotoBins)
    real*8,value :: Tgas
    real*8::tau,n(nspec)
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
       tau = 0d0
       !loop on species
       do i=1,nPhotoRea
          !calc opacity as column_density * cross_section
          idx = photoPartners(i)
          tau = tau + num2col(x(idx),n(:)) * photoBinJTab(i,j)
       end do
       krome_get_opacity(j) = tau !store
    end do

  end function krome_get_opacity

  !*****************************
  !get the opacity tau corresponding to the x(:)
  ! chemical composition. The column density
  ! is computed using the size of the cell (csize)
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity_size(x,Tgas,csize)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_dust
    implicit none
    real*8 :: x(nmols),krome_get_opacity_size(nPhotoBins)
    real*8,value :: Tgas,csize
    real*8::n(nspec),energy,tau
    integer::i,j,idx

    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    !loop on frequency bins
    do j=1,nPhotoBins
       tau = 0d0
       !loop on species
       do i=1,nPhotoRea
          !calc opacity as column_density * cross_section
          !where column_density is density*cell_size
          idx = photoPartners(i)
          tau = tau + x(idx) * photoBinJTab(i,j)
       end do

#IFKROME_dust_opacity
       energy = photoBinEmid(j)
       do i=1,ndust
          tau = tau + pi*krome_dust_asize2(i)*xdust(i) * get_Qabs(energy,i)
       end do
#ENDIFKROME_dust_opacity
       krome_get_opacity_size(j) = tau * csize !store
    end do

#IFKROME_useH2pd
    krome_get_opacity_size(:) = krome_get_opacity_size(:) &
         + kpd_H2_xsec(Tgas) * n(idx_H2) * csize
#ENDIFKROME_useH2pd

  end function krome_get_opacity_size

  !*****************************
  !get the opacity tau corresponding to the x(:)
  ! chemical composition. The column density
  ! is computed using the size of the cell (csize).
  ! Dust is included using dust-to-gas mass ratio (d2g).
  ! You should load the dust tables with the subroutine
  ! krome_load_dust_opacity(fileName).
  ! An array of size krome_nPhotoBins is returned.
  function krome_get_opacity_size_d2g(x,Tgas,csize,d2g)
    use krome_commons
    use krome_constants
    use krome_photo
    use krome_subs
    use krome_dust
    use krome_getphys
    implicit none
    real*8 :: x(nmols),krome_get_opacity_size_d2g(nPhotoBins)
    real*8,value :: Tgas,csize,d2g
    real*8::n(nspec),energy,tau,m(nspec),mgas
    integer::i,j,idx

    m(:) = get_mass()
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    mgas = sum(n(1:nmols)*m(1:nmols))

    !loop on frequency bins
    do j=1,nPhotoBins
       tau = 0d0
       !loop on species
       do i=1,nPhotoRea
          !calc opacity as column_density * cross_section
          !where column_density is density*cell_size
          idx = photoPartners(i)
          tau = tau + x(idx) * photoBinJTab(i,j)
       end do

       !sum dust opacity from interpolated table
       tau = tau + d2g*mgas*opacityDust(j)

       krome_get_opacity_size_d2g(j) = tau * csize !store
    end do

  end function krome_get_opacity_size_d2g

  !*********************
  !scale radiation intensity with opacity assuming a given
  ! cell size and gas composition
!  subroutine krome_opacity_scale_size(csize,n,Tgas)
!    use krome_commons
!    implicit none
!    real*8::csize,n(nmols),xscale(nPhotoBins),Tgas
!
!    xscale(:) = krome_get_opacity_size(n(:),Tgas,csize)
!    xscale(:) = exp(-xscale(:))
!    call krome_photoBin_scale_array(xscale(:))
!
!  end subroutine krome_opacity_scale_size

  !*********************
  !scale radiation intensity with opacity assuming a given
  ! cell size and gas composition
  subroutine krome_opacity_scale_size(csize,n,Tgas) #KROME_bindC
    use krome_commons
    implicit none
#IFKROME_useBindC
    real(kind=c_double),value::csize,Tgas
    real(kind=c_double) :: n(nmols)
#ELSEKROME_useBindC
    real*8::csize,n(nmols),Tgas
#ENDIFKROME_useBindC
    real*8 :: xscale(nPhotoBins)

    xscale(:) = krome_get_opacity_size(n(:),Tgas,csize)
    xscale(:) = exp(-xscale(:))
    call krome_photoBin_scale_array(xscale(:))
  end subroutine krome_opacity_scale_size


  !*******************************
  !load a frequency-dependent opacity table stored in fname file,
  ! column 1 is energy or wavelenght in un units of unitEnergy
  ! (default eV), column 2 is opacity in cm2/g.
  ! opacity is interpolated over the current photo-binning.
#IFKROME_useBindC
  subroutine krome_load_opacity_table() bind(C)
    !C-Fortran interoperability has problems dealing with strings
    !Assume the opacity file name is 'opacityDust.dat'
#ELSEKROME_useBindC
  subroutine krome_load_opacity_table(fname, unitEnergy)
#ENDIFKROME_useBindC
    use krome_commons
    use krome_constants
    use krome_photo
    implicit none
#IFKROME_useBindC
    character(len=80)::fname
#ELSEKROME_useBindC    
    character(len=*)::fname
    character(len=*),optional::unitEnergy
#ENDIFKROME_useBindC
    character*10::eunit
    
    !read energy unit optional argument
    eunit = "eV" !default is eV
#IFKROME_useBindC
    fname = "opacityDust.dat"
#ELSEKROME_useBindC    
    if(present(unitEnergy)) then
      eunit = trim(unitEnergy)
    end if
#ENDIFKROME_useBindC    

    call load_opacity_table(fname, eunit)

end subroutine krome_load_opacity_table

   !*******************************
  !load a frequency-dependent opacity table stored in fname file,
  ! column 1 is energy or wavelenght in un units of unitEnergy
  ! (default eV), column 2 is opacity in cm2/g.
  ! opacity is interpolated over the current photo-binning.
  subroutine krome_load_opacity_table(fname, unitEnergy)
    use krome_commons
    use krome_constants
    use krome_photo
    implicit none
    integer,parameter::ntmp=int(1e5)
    character(len=*)::fname
    character(len=*),optional::unitEnergy
    character*10::eunit
    integer::ios,icount,iR,iL,i,j,fileUnit
    real*8::wl,opac,fL,fR,kk,dE
    real*8::wls(ntmp),opacs(ntmp)
    real*8,allocatable::energy(:),kappa(:)

    !read energy unit optional argument
    eunit = "eV" !default is eV
    if(present(unitEnergy)) then
      eunit = trim(unitEnergy)
    end if

    call load_opacity_table(fname, eunit)

end subroutine krome_load_opacity_table

  ! ******************************
  ! load absorption data data from file, cm2/g
  subroutine krome_load_average_kabs()
    use krome_photo
    implicit none

    call find_Av_load_kabs()

  end subroutine krome_load_average_kabs

  ! *******************************
  ! use linear least squares and the current Jflux distribution
  ! to return G0 and Av.
  ! x(:) are the abundances (use for mean molecular weight)
  ! and d2g is dust to gas mass ratio
  subroutine krome_find_G0_Av(G0, Av, x, d2g)
    use krome_commons
    use krome_photo
    implicit none
    real*8,intent(out)::G0, Av
    real*8,intent(in)::d2g, x(nmols)
    real*8::n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0

    call estimate_G0_Av(G0, Av, n(:), d2g)

  end subroutine krome_find_G0_Av

  !*******************************
  !dump the Jflux profile to the file
  ! with unit number nfile
  subroutine krome_dump_Jflux(nfile) #KROME_bindC
    use krome_commons
    implicit none
    integer::i
#KROME_integer_value :: nfile

    do i=1,nPhotoBins
       write(nfile,*) photoBinEmid(i),photoBinJ(i)
    end do

  end subroutine krome_dump_Jflux

  !**********************
  !set the velocity for line broadening, cm/s
  subroutine krome_set_lineBroadeningVturb(vturb)
    use krome_constants
    use krome_commons
    implicit none
    real*8::vturb

    broadeningVturb2 = vturb**2

  end subroutine krome_set_lineBroadeningVturb

#ENDIFKROME

  !***************************
  !alias for coe in krome_subs
  ! returns the coefficient array of size krome_nrea
  ! for a given Tgas
  function krome_get_coef(Tgas,x)
    use krome_commons
    use krome_subs
    use krome_tabs
    real*8 :: krome_get_coef(nrea),x(nmols)
    real*8,value:: Tgas
    real*8::n(nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

#IFKROME_hasStoreOnceRates
    call makeStoreOnceRates(n(:))
#ENDIFKROME

    krome_get_coef(:) = coe(n(:))

  end function krome_get_coef

  !****************************
  !get the mean molecular weight from
  ! mass fractions
  function krome_get_mu_x(xin) #KROME_bindC
    use krome_commons
    implicit none
#KROME_double :: xin(nmols), krome_get_mu_x
    real*8::n(nmols)
    n(:) = krome_x2n(xin(:),1d0)
    krome_get_mu_x = krome_get_mu(n(:))
  end function krome_get_mu_x

  !****************************
  !return the adiabatic index from mass fractions
  ! and temperature in K
  function krome_get_gamma_x(xin,inTgas) #KROME_bindC
    use krome_commons
    implicit none
#KROME_double_value :: inTgas
#KROME_double :: xin(nmols), krome_get_gamma_x
    real*8::x(nmols),Tgas,rhogas

    Tgas = inTgas
    x(:) = krome_x2n(xin(:),1d0)
    krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

  end function krome_get_gamma_x


  !***************************
  !normalize mass fractions and
  ! set charge to zero
  subroutine krome_consistent_x(x) #KROME_bindC
    use krome_commons
    use krome_constants
    implicit none
#KROME_double :: x(nmols)
    real*8::isumx,sumx,xerr,imass(nmols),ee

    !1. charge consistency
    imass(:) = krome_get_imass()

#KROME_zero_electrons

    ee = sum(krome_get_charges()*x(:)*imass(:))
    ee = max(ee*e_mass,0d0)
#KROME_electrons_balance

    !2. mass fraction consistency
    sumx = sum(x)

    !NOTE: uncomment here if you want some additional control
    !conservation error threshold: rise an error if above xerr
    !xerr = 1d-2
    !if(abs(sum-1d0)>xerr) then
    !   print *,"ERROR: some problem with conservation!"
    !   print *,"|sum(x)-1|=",abs(sum-1d0)
    !   stop
    !end if

    isumx = 1d0/sumx
    x(:) = x(:) * isumx

  end subroutine krome_consistent_x

  !*********************
  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    use krome_commons
    implicit none
    real*8 :: n(nmols),krome_n2x(nmols)
    real*8,value :: rhogas

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    use krome_commons
    implicit none
    real*8 :: x(nmols),krome_x2n(nmols)
    real*8,value :: rhogas

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !******************
  !returns free-fall time using the number density
  ! abundances of array x(:)
  function krome_get_free_fall_time(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time
    real*8::x(:),n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_free_fall_time = get_free_fall_time(n(:))

  end function krome_get_free_fall_time

  !******************
  !returns free-fall time using the total mass density
  !  of gas, rhogas (g/cm3)
  function krome_get_free_fall_time_rho(rhogas)
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time_rho
    real*8::rhogas

    krome_get_free_fall_time_rho = get_free_fall_time_rho(rhogas)

  end function krome_get_free_fall_time_rho

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt) #KROME_bindC
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    use krome_gadiab
    implicit none
#KROME_double :: x(nmols), Tgas
#KROME_double_value :: dt
    real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

#IFKROME_use_thermo
    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
         * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas
#ENDIFKROME

  end subroutine krome_thermo

#IFKROME_use_heating
  !*************************
  !get heating (erg/cm3/s) for a given species
  ! array x(:) and Tgas
  function krome_get_heating(x,inTgas) #KROME_bindC
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
#KROME_double_value :: inTgas
#KROME_double :: x(nmols), krome_get_heating
    real*8::Tgas,k(nrea),nH2dust,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    k(:) = coe(n(:))
    nH2dust = 0d0
    krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
  end function krome_get_heating

  !*****************************
  ! get an array containing individual heatings (erg/cm3/s)
  ! the array has size krome_nheats. see heatcool.gps
  ! for index list
  function krome_get_heating_array(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k(nrea),nH2dust
    real*8 :: x(nmols),krome_get_heating_array(nheats)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    k(:) = coe(n(:))
    Tgas = inTgas
    nH2dust = 0d0
    krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

  end function krome_get_heating_array

#ENDIFKROME

#IFKROME_use_cooling
  !*************************
  !get cooling (erg/cm3/s) for x(:) species array
  ! and Tgas
  function krome_get_cooling(x,inTgas) #KROME_bindC
    use krome_cooling
    use krome_commons
    implicit none
#KROME_double_value :: inTgas
#KROME_double :: x(nmols), krome_get_cooling
    real*8::Tgas,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    krome_get_cooling = cooling(n,Tgas)
  end function krome_get_cooling

  !*****************************
  ! get an array containing individual coolings (erg/cm3/s)
  ! the array has size krome_ncools. see heatcool.gps
  ! for index list
  function krome_get_cooling_array(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas
    real*8 :: x(nmols),krome_get_cooling_array(ncools)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    Tgas = inTgas
    krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

  end function krome_get_cooling_array

  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n) #KROME_bindC
    use krome_cooling
    implicit none
#KROME_double :: n(krome_nmols)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in) #KROME_bindC
    use krome_cooling
    use krome_commons
    implicit none
#KROME_double :: n(nmols)
#KROME_double_value :: Tgas
    real*8::x(nspec)
#IFKROME_useBindC
#KROME_integer_value :: nfile_in
#ELSEKROME_useBindC
#KROME_integer, optional :: nfile_in
#ENDIFKROME_useBindC
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
#IFKROME_useBindC
#ELSEKROME_useBindC
    if(present(nfile_in)) nfile = nfile_in
#ENDIFKROME_useBindC
    call dump_cool(x(:),Tgas,nfile)

  end subroutine krome_dump_cooling

#ENDIFKROME

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  subroutine krome_conserveLin_x(x,ref) #KROME_bindC
    use krome_commons
    use krome_subs
    implicit none
#KROME_double :: x(nmols),ref(natoms)

    call conserveLin_x(x(:),ref(:))

  end subroutine krome_conserveLin_x

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  function krome_conserveLinGetRef_x(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

    krome_conserveLinGetRef_x(:) = &
         conserveLinGetRef_x(x(:))

  end function krome_conserveLinGetRef_x

  !*************************
  !force conservation to array x(:)
  !using xi(:) as initial abundances.
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
    real*8::n(krome_nspec),ni(krome_nspec)

    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !get the adiabatic index for x(:) species abundances
  ! and Tgas.
  ! alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas) #KROME_bindC
    use krome_subs
    use krome_commons
    use krome_gadiab
#KROME_double_value :: Tgas
#KROME_double :: x(nmols), krome_get_gamma
    real*8::n(nspec)
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies.
  ! alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    integer :: krome_get_zatoms(nmols)
    integer::zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)

  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from
  ! number density and mass density.
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
#KROME_double :: x(nmols), krome_get_mu
    real*8::n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    krome_get_mu = get_mu(n(:))
  end function krome_get_mu

  !***************************
  !get the names of the reactions as a
  ! character*50 array of krome_nrea
  ! elements
  !!#KROME_bindC !! cannot yet be called from C
  function krome_get_rnames()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    character*50 :: krome_get_rnames(nrea)

    krome_get_rnames(:) = get_rnames()

  end function krome_get_rnames

  !*****************
  !get an array of double containing the masses in g
  ! of the species.
  ! alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(nspec)
#KROME_double :: krome_get_Hnuclei, x(nmols)
    n(:) = 0d0
    n(1:nmols) = x(:)

    krome_get_Hnuclei = get_Hnuclei(n(:))

  end function krome_get_Hnuclei

  !*****************
  !get an array of size krome_nmols containing the
  ! charges of the species.
  ! alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 and size krome_nmols
  ! containing the names of all the species.
  ! alias for get_names
  !! #KROME_bindC !! cannot yet be called from C
  function krome_get_names()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    character*16 :: krome_get_names(nmols)
    character*16::tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !********************
  !get space-separated header of chemical species
  function krome_get_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    #KROME_names_header_define
    character*16::tmp(nspec)
    integer::i

    tmp(:) = get_names()

    krome_get_names_header = ""
    do i=1,nmols
       krome_get_names_header = trim(krome_get_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_names_header

  !********************
  !get space-separated header of coolings
  function krome_get_cooling_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    #KROME_cooling_names_header_define
    character*16::tmp(ncools)
    integer::i

    tmp(:) = get_cooling_names()

    krome_get_cooling_names_header = ""
    do i=1,ncools
       if(trim(tmp(i))=="") cycle
       krome_get_cooling_names_header = trim(krome_get_cooling_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_cooling_names_header

  !********************
  !get space-separated header of heatings
  function krome_get_heating_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    #KROME_heating_names_header_define
    character*16::tmp(nheats)
    integer::i

    tmp(:) = get_heating_names()

    krome_get_heating_names_header = ""
    do i=1,nheats
       if(trim(tmp(i))=="") cycle
       krome_get_heating_names_header = trim(krome_get_heating_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_heating_names_header

  !*****************
  !get the index of the species with name name.
  ! alias for get_index
  !!#KROME_bindC !! cannot yet be called from C
  function krome_get_index(name)
    use krome_subs
    implicit none
#KROME_integer :: krome_get_index
    character*(*) :: name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n) #KROME_bindC
    use krome_commons
#KROME_double :: krome_get_rho, n(nmols)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho

  !*************************
  !scale the abundances of the metals contained in n(:)
  ! to Z according to Asplund+2009.
  ! note that this applies only to neutral atoms.
  subroutine krome_scale_Z(x,Z) #KROME_bindC
    use krome_commons
    use krome_getphys
#KROME_double :: x(nmols)
#KROME_double_value :: Z
    real*8::Htot,n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0

#KROME_scaleZ

  end subroutine krome_scale_Z

  !*************************
  !set the total metallicity
  ! in terms of Z/Z_solar
  subroutine krome_set_Z(xarg) #KROME_bindC
    use krome_commons
#KROME_double_value :: xarg

    total_Z = xarg

  end subroutine krome_set_Z

  !*************************
  !set D is in terms of D_solar (D/D_sol).
  subroutine krome_set_dust_to_gas(xarg) #KROME_bindC
     use krome_commons
#KROME_double_value :: xarg

   dust2gas_ratio = xarg

  end subroutine

  !*************************
  !set the clumping factor
  subroutine krome_set_clump(xarg) #KROME_bindC
    use krome_commons
#KROME_double_value :: xarg

    clump_factor = xarg

  end subroutine krome_set_clump


  !***********************
  !get the number of electrons assuming
  ! total neutral charge (cations-anions)
  function krome_get_electrons(x) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_getphys
#KROME_double :: x(nmols), krome_get_electrons
    real*8::n(nspec)
    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_electrons = get_electrons(n(:))
  end function krome_get_electrons

  !**********************
  !print on screen the first nbest highest reaction fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest) #KROME_bindC
    use krome_subs
    use krome_commons
    implicit none
#KROME_double :: xin(nmols)
#KROME_double_value :: Tgas
    real*8::x(nmols),n(nspec)
#KROME_integer_value :: nbest
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)

  end subroutine krome_print_best_flux

  !*********************
  !print only the highest fluxes greater than a fraction frac
  ! of the maximum flux
  subroutine krome_print_best_flux_frac(xin,Tgas,frac) #KROME_bindC
    use krome_subs
    use krome_commons
    implicit none
#KROME_double :: xin(nmols)
#KROME_double_value :: Tgas,frac
    real*8::n(nspec)
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_frac(n,Tgas,frac)

  end subroutine krome_print_best_flux_frac

  !**********************
  !print the highest nbest fluxes for reactions involving
  !a given species using the index idx_find (e.g. krome_idx_H2)
  subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find) #KROME_bindC
    use krome_subs
    use krome_commons
    implicit none
#KROME_double :: xin(nmols)
#KROME_double_value :: Tgas
    real*8::n(nspec)
#KROME_integer_value :: nbest,idx_find
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_spec(n,Tgas,nbest,idx_find)
  end subroutine krome_print_best_flux_spec

  !*******************************
  !get an array of size krome_nrea with
  ! the fluxes of all the reactions in cm-3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8 :: krome_get_flux(nrea),n(nmols)
    real*8,value :: Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*****************************
  !store the fluxes to the file unit ifile
  ! using the chemical composition x(:), and the
  ! gas temperature Tgas. xvar is th value of an
  ! user-defined independent variable that
  ! can be employed for plots.
  ! the file columns are as follow
  ! rate number, xvar, absolute flux,
  !  flux/maxflux, flux fraction wrt total,
  !  reaction name (*50 string)
  subroutine krome_explore_flux(x,Tgas,ifile,xvar) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
#KROME_double :: x(nmols)
#KROME_double_value :: Tgas,xvar
    real*8::flux(nrea),fluxmax,sumflux,n(nspec)
#KROME_integer_value :: ifile
    integer::i
    character*50::rname(nrea)

    !get reaction names
    rname(:) = get_rnames()
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    !get fluxes
    flux(:) = get_flux(n(:), Tgas)
    fluxmax = maxval(flux) !maximum flux
    sumflux = sum(flux) !sum of all the fluxes
    !loop on reactions
    do i=1,nrea
       write(ifile,'(I8,5E17.8e3,a3,a50)') i,xvar,Tgas,flux(i),&
            flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
    end do
    write(ifile,*)

  end subroutine krome_explore_flux

  !*********************
  !get nulcear qeff for the reactions
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: krome_get_qeff(nrea)

    krome_get_qeff(:) = get_qeff()

  end function krome_get_qeff

#IFKROME_useStars

  !**************************
  !alias for stars_coe in krome_star
  function krome_stars_coe(x,rho,Tgas,y_in,zz_in)
    use krome_commons
    use krome_stars
    use krome_subs
    use krome_getphys
    implicit none
    real*8,value :: rho,Tgas
    real*8 :: x(nmols)
    real*8 :: krome_stars_coe(nrea)
    real*8::n(nspec)
    integer::ny
    real*8,allocatable::y(:)
    integer,allocatable::zz(:)
    real*8, optional::y_in(:)
    integer, optional::zz_in(:)

    !check if extened abundances and zatom array are present
    if(present(y_in)) then
       ny = size(y_in)
       allocate(y(ny), zz(ny))
       y(:) = y_in(:)
       zz(:) = zz_in(:)
    else
       allocate(y(nmols),zz(nmols))
       zz(:) = get_zatoms()
       y(:) = x(:)
    end if

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_stars_coe(:) = stars_coe(n(:),rho,Tgas,y(:),zz(:))
    deallocate(y,zz)
  end function krome_stars_coe

  !********************************
  !alias for stars_energy
  function krome_stars_energy(x,rho,Tgas,k)
    use krome_commons
    use krome_stars
    implicit none
    real*8 :: x(nmols),k(nrea)
    real*8,value :: rho,Tgas
    real*8 :: krome_stars_energy(nrea)
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_stars_energy(:) = stars_energy(n(:),rho,Tgas,k(:))

  end function krome_stars_energy

#ENDIFKROME

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile) #KROME_bindC
    use krome_commons
#KROME_double :: n(nmols)
#KROME_double_value :: Tgas
    real*8::flux(nrea)
#KROME_integer_value :: nfile
    integer::i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
       write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  !dump all the evaluation of the coefficient rates in
  ! the file funit, in the range inTmin, inTmax, using
  ! imax points
  subroutine krome_dump_rates(inTmin,inTmax,imax,funit) #KROME_bindC
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
#KROME_integer_value :: funit,imax
#KROME_double_value :: inTmin,inTmax
    real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

    Tmin = log10(inTmin)
    Tmax = log10(inTmax)

    n(:) = 1d-40
    do i=1,imax
       Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
       n(idx_Tgas) = Tgas
       k(:) = coe(n(:))
       do j=1,nrea
          write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
       end do
       write(funit,*)
    end do


  end subroutine krome_dump_rates

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas) #KROME_bindC
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer::i,charges(nspec)
#KROME_double :: x(nmols)
#KROME_double_value :: Tgas
    real*8::masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
       print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info

  !*****************************
  subroutine krome_set_mpi_rank(xarg) #KROME_bindC
    use krome_commons
    implicit none
    #KROME_integer_value :: xarg
    krome_mpi_rank=xarg
  end subroutine krome_set_mpi_rank

#IFKROME_usePhotoBins
#IFKROME_useBindC
  !*****************************
  function krome_get_photoIntensity(energy) bind(C)
    use krome_commons
    use krome_photo
    implicit none
    real(kind=c_double),value::energy
    real(kind=c_double),target::krome_get_photoIntensity
    krome_get_photoIntensity = get_photoIntensity(energy)
  end function krome_get_photoIntensity
#ENDIFKROME_useBindC
#ENDIFKROME


  !**************************
  function krome_get_jacobian(j,x,Tgas)
    use krome_ode
    use krome_commons
    implicit none
    integer, value :: j
    real*8,value :: Tgas
    real*8 :: x(nmols),krome_get_jacobian(nspec)
    integer::ian, jan, i
    real*8::tt, n(nspec)
    real*8::pdj(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = tgas

    tt = 0d0
    ian = 0
    jan = 0

    call jes(nspec, tt, n, j, ian, jan, pdj)
    krome_get_jacobian(:) = pdj(:)

  end function krome_get_jacobian

end module krome_user
