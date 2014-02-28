subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
  character(LEN=80)::filename
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  

  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,icpu,mcpu,info
  real(dp)::d,u,v,w,A,B,C,e
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  character(LEN=80)::id='$Id: output_hydro.f90,v 1.9 2013/07/20 17:57:45 troels_h Exp $'

  call print_id(id)

  if(verbose)write(*,*)'Entering backup_hydro'

  !-----------------------------------
  ! Output hydro vars in file. Max 128 threads output simultaneously.
  !-----------------------------------  
  mcpu = (ncpu_dump-1)/128+1                       ! number of times we have to dump
  do icpu=0,mcpu-1                            ! Nr of times we dump
     if (myid==1 .and. mcpu > 3) print '(a,i3,a,i3)', 'Dumping MHD batch ',icpu+1,' out of ',mcpu
     if (modulo(myid-1,mcpu)==icpu .and. ngrid_current>0) then
        ilun=10
        call title(myid,nchar)
        fileloc=TRIM(filename)//TRIM(nchar)
        open(unit=ilun,file=fileloc,form='unformatted')
        write(ilun)ncpu_dump
        write(ilun)nvar+3
        write(ilun)ndim
        write(ilun)nlevelmax
        write(ilun)nboundary
        write(ilun)gamma
        do ilevel=1,nlevelmax
           do ibound=1,nboundary+ncpu_dump
              if(ibound<=ncpu_dump)then
                 ncache=numbl(ibound,ilevel)
                 istart=headl(ibound,ilevel)
              else
                 ncache=numbb(ibound-ncpu_dump,ilevel)
                 istart=headb(ibound-ncpu_dump,ilevel)
              end if
              write(ilun)ilevel
              write(ilun)ncache
              if(ncache>0)then
                 allocate(xdp(1:ncache))
                 ! Loop over cells
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    do ivar=1,4
                       igrid=istart
                       if(ivar==1)then ! Write density
                          do i=1,ncache
                             xdp(i)=uold(igrid+iskip,1)
                             igrid=next(igrid)
                          end do
                       else ! Write velocity field
                          if (write_conservative) then
                             do i=1,ncache
                                xdp(i)=uold(igrid+iskip,ivar)
                                igrid=next(igrid)
                             end do
                          else
                             do i=1,ncache
                                xdp(i)=uold(igrid+iskip,ivar)/uold(igrid+iskip,1)
                                igrid=next(igrid)
                             end do
                          endif
                       endif
                       write(ilun)xdp
                    end do
                    do ivar=6,8 ! Write left B field
                       igrid=istart
                       do i=1,ncache
                          xdp(i)=uold(igrid+iskip,ivar)
                          igrid=next(igrid)
                       end do
                       write(ilun)xdp
                    end do
                    do ivar=nvar+1,nvar+3 ! Write right B field
                       igrid=istart
                       do i=1,ncache
                          xdp(i)=uold(igrid+iskip,ivar)
                          igrid=next(igrid)
                       end do
                       write(ilun)xdp
                    end do
                    igrid=istart
                    do i=1,ncache ! Write pressure
                       if (write_conservative) then
                          xdp(i)=uold(igrid+iskip,5)
                       else
                          d=uold(igrid+iskip,1)
                          u=uold(igrid+iskip,2)/d
                          v=uold(igrid+iskip,3)/d
                          w=uold(igrid+iskip,4)/d
                          A=0.5*(uold(igrid+iskip,6)+uold(igrid+iskip,nvar+1))
                          B=0.5*(uold(igrid+iskip,7)+uold(igrid+iskip,nvar+2))
                          C=0.5*(uold(igrid+iskip,8)+uold(igrid+iskip,nvar+3))
                          e=uold(igrid+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
                          xdp(i)=(gamma-1d0)*e
                       endif
                       igrid=next(igrid)
                    end do
                    write(ilun)xdp
                    do ivar=9,nvar ! Write passive scalars if any
                       igrid=istart
                       do i=1,ncache
                          xdp(i)=uold(igrid+iskip,ivar)!/uold(igrid+iskip,1)
                          igrid=next(igrid)
                       end do
                       write(ilun)xdp
                    end do
                 end do
                 deallocate(xdp)
              end if
           end do
        end do
        close(ilun)
     end if
#ifndef WITHOUTMPI
     call MPI_BARRIER(mpi_comm_use,info)
#endif
  end do
     
end subroutine backup_hydro

subroutine dump_hydro(var,filename)
  use amr_commons
  use hydro_commons
  implicit none
  real(kind=dp) var(:)
  character(LEN=*)::filename
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  

  integer::i,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,icpu,mcpu,info,nvar1
  real(dp)::d,u,v,w,A,B,C,e
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar,nmyid
  character(LEN=80)::fileloc,filedir

  if(verbose)write(*,*)'Entering dump_hydro'

  call title(ifout,nchar)
  filedir=TRIM(datadir)//'/'//'output_'//TRIM(nchar)//'/'

  !-----------------------------------
  ! Output hydro vars in file. Max 128 threads output simultaneously.
  !-----------------------------------  
  mcpu = (ncpu_dump-1)/128+1                       ! number of times we have to dump
  nvar1 = 1
  do icpu=0,mcpu-1                            ! Nr of times we dump
     if (modulo(myid-1,mcpu)==icpu) then
        ilun=10
        call title(ifout,nchar)
        call title(myid,nmyid)
        fileloc=TRIM(filedir)//TRIM(filename)//TRIM(nchar)//'.out'//TRIM(nmyid)
        print*,myid,TRIM(fileloc)
        open(unit=ilun,file=fileloc,form='unformatted')
        write(ilun)ncpu_dump
        write(ilun)nvar1
        write(ilun)ndim
        write(ilun)nlevelmax
        write(ilun)nboundary
        write(ilun)gamma
        do ilevel=1,nlevelmax
           do ibound=1,nboundary+ncpu_dump
              if(ibound<=ncpu_dump)then
                 ncache=numbl(ibound,ilevel)
                 istart=headl(ibound,ilevel)
              else
                 ncache=numbb(ibound-ncpu_dump,ilevel)
                 istart=headb(ibound-ncpu_dump,ilevel)
              end if
              write(ilun)ilevel
              write(ilun)ncache
              if(ncache>0)then
                 allocate(xdp(1:ncache))
                 ! Loop over cells
                 do ind=1,twotondim
                    iskip=ncoarse+(ind-1)*ngridmax
                    igrid=istart
                    do i=1,ncache
                       xdp(i)=i
                       !xdp(i)=var(igrid+iskip)
                       igrid=next(igrid)
                    end do
                    write(ilun)xdp
                 end do
                 deallocate(xdp)
              end if
           end do
        end do
        close(ilun)
     end if
#ifndef WITHOUTMPI
     call MPI_BARRIER(mpi_comm_use,info)
#endif
  end do

end subroutine dump_hydro
