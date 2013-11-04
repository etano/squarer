      function vsum0(m,i,p,t,n)
      include 'mach.p'
      real*8 vsum0,p,t
      integer m,i,n,j
      dimension i(m),p(m),t(n,4)
      vsum0=0.
      do j=1,m
        vsum0 = vsum0 + t(i(j),1) 
     &                + p(j)*(t(i(j),2)
     &                        + p(j)*(t(i(j),3)
     &                                + p(j)*t(i(j),4)))
      end do
      return
      end
