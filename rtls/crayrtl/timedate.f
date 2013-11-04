      subroutine timedate(i)
      character i*26,i1*8,i2*8,name*20
      call clock(i1)
      i(1:8)=i1(1:8)    
      call date(i2)
      i(9:9)=' '
      i(10:17)=i2(1:8)  
      call gethost(name)
      i(18:26)=name(1:9)
      write (6,*)' time = ',i(1:17) ,' host=',name
      return
      end
