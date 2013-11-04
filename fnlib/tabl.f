      subroutine tabl(m,i,p,t,n,result)
      include 'mach.p'
      integer m,i(m),n,j
      real*8 p(m),t(n,4),result(m)
      do j=1,m
        result(j)=t(i(j),1) 
     &                + p(j)*(t(i(j),2)
     &                        + p(j)*(t(i(j),3)
     &                                + p(j)*t(i(j),4)))
      end do
      return
      end
