      subroutine sumavv(iwrt)
      include 'mach.p'

      integer iwrt,ify,i
      include 'caver.cm'
      real*8 avcu0,avbl(maver),avcu,error

c increments the averges counters. iwrt nonzero means no prints
      if(nblock.lt.0) then
      if(iwrt.eq.0) write (*,*)'nblock not initialized', nblock
      stop
      endif
      if(nblock.eq.0) then
            write (51,'(a7,100(1H",a16,2H" ))')'# '
     .        ,(qaname(i),i=1,naver)
      endif


      ify=0 ! check if any work done this block
      do i=1,naver
        if(anorm(i).ne.0.d0) ify=1
      enddo
      if(ify.eq.0) return

c write out data on unit 62
!      write (62,62) (avtemp(i),i=1,maver)
!      write (62,62) (anorm(i),i=1,maver)
!62     format(6e18.10)
!      write (62,63) (qaname(i),i=1,maver)
!      write (62,63) (aunits(i),i=1,maver)
!63     format(6a8)

      nblock=nblock+1
      write (62) naver,nblock,(avtemp(i),anorm(i),i=1,naver)

      if(iwrt.eq.0) write (6,33) nblock
33    format(/,1x,20(1h*),14hfinished block,i5,34(1h*)//
     . ,5x,4hname,5x,5hunits,1x,2(5hblock,7x),3(5hcuml.,7x)
     .  /18x,2(7haverage,5x,4hnorm,8x),5herror)
      do i=1,naver
      if(anorm(i).ne.0.0d0) then

         ablock(i)=ablock(i)+1.0d0
         if(ansum(i).ne.0.0d0) avcu0=avsum(i)/ansum(i)
         ansum(i)=ansum(i)+anorm(i)
         avsum(i)=avsum(i)+avtemp(i)
         avbl(i)=avtemp(i)/anorm(i)

         if(ansum(i).ne.0.0) then
             avcu=avsum(i)/ansum(i)
         else
           avcu=0.0d0
         endif

         if(ablock(i).gt.1.9.and.ansum(i).ne.0.d0) then
            avsq(i)=avsq(i)+(1.0d0-anorm(i)/ansum(i))*
     .                   (avbl(i)-avcu0)**2*anorm(i)
            error=sqrt(avsq(i)/(ansum(i)*(ablock(i)-1.d0)))
         else
             error=0.0d0
         endif

       if(iwrt.eq.0)
     .    write (6,30) i,qaname(i),aunits(i),avbl(i),anorm(i)
     .                 ,avcu,ansum(i),error
30    format(1h ,i2,1x,a9,1x,a3,5e12.5)
      else
      avbl(i)=0.d0
      endif

      anorm(i)=0.0d0
      avtemp(i)=0.0d0
      enddo
      write (51,'(100e14.6)')(avbl(i),i=1,naver)

      if(iwrt.eq.0) write (6,'(/)')
      call flush(62)
      return
      end 
