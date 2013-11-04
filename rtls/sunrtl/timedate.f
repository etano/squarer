      subroutine timedate(i)
      character fdate*24,name*12,i*26
      call hostnm(name)
c delete the day of week field?
      i(1:19)=fdate()
      i(20:26)=name(1:7)
      write (6,*)' time = ',i,' host ',name
      return
      end
