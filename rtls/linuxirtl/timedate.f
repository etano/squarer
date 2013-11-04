      subroutine timedate(i)
      character fdate*24,name*12,i*26
      i= '                   '         
      name='intel   '
!     call hostnm(name)
c delete the day of week field?
!     i(1:19)=fdate()
      i(18:19)='  '
      i(20:26)=name(1:7)
!     call date(i(1:9))
!     call time(i(10:17))
!     write (6,'(a9,1x,a8,1x,a7)')i(1:9),i(10:17),name(1:7)
      call date_and_time(i(1:8),i(9:18))
      write (6,'(a4,1x,a2,1x,a2,1x,a2,1x,2a,1x,6a,1x,a7)')
     & i(1:4),i(5:6),i(7:8),i(9:10),i(11:12),i(13:18),name(1:7)
      end
