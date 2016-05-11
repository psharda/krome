      program test
      USE frt_cf3_mod

      real :: GHab
      real(kind=8) :: Tem, Den, Plw, Ph1, Pg1, Pc6, cfun1, hfun1, cfun0, hfun0

      mode = 0
      call frtInitCF(mode,'cf_table.I2.dat')
      if(mode .ne. 0) then
         write(0,*) 'Error in reading data file: ', mode
         stop
      endif

      !Den = 1.0e-3
      !Plw = 2.11814e-13
      !Ph1 = 1.08928e-13
      !Pg1 = 2.76947e-14
      !Pc6 = 1.03070e-17
      Den = 1e4
      GHab=1.7
      Plw = GHab*4e-11
      Ph1 = 1.08928e-13
      Pg1 = 2.76947e-14
      Pc6 = 1.03070e-17

      alt = 1.0
      do while(alt .le. 9.0)
         Tem = 10.0**alt

         call frtGetCF(Tem,Den,1.0_8,Plw,Ph1,Pg1,Pc6,cfun0,hfun0,ierr)
         if(ierr .ne. 0) then
            write(0,*) 'Error in table call: ', ierr
            stop
         endif

         call frtGetCF(Tem,Den,1.0_8,Plw,Ph1,Pg1,Pc6,cfun1,hfun1,ierr)
         if(ierr .ne. 0) then
            write(0,*) 'Error in table call: ', ierr
            stop
         endif

         write(6,9) alt, cfun0, hfun0, cfun1, hfun1
 9       format(F3.1,1P,4(1X,E10.3))

         alt = alt + 0.1

      enddo

      end

