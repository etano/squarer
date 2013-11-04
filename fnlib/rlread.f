      function rlread(p)
      include 'mach.p'
      integer ln,ib,i
      real*8 rlread
      character p*(*),blank*28,both*28
      data blank/'                    '/ 
      ln=index(p,' ')-1
!check for a decimal place
      ib=0
      do i=1,ln
       if(p(i:i).eq.'.') ib=1
      enddo
      if(ib.eq.1) then
      both=blank(1:24-ln)//p(1:ln)
      else
       both=blank(1:23-ln)//p(1:ln)//'.'
      endif
      read(both,1) rlread
1     format(e24.14)
      if(ib.eq.0) write (*,*)' . appended in rlread ',rlread
      return
      end
