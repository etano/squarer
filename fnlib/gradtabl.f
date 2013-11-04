      subroutine gradtabl(m,i,p,t,fp,n,dxdr,r,pact)
      include 'mach.p'
c computes  first derivative of spline only (1/r*du/dr)
      integer m,i(m),n,k
      real*8 p(m),t(n,4),fp(m),r(m),dxdr(m),pact
      pact=0.
      do k=1,m
       fp(k) = dxdr(k)*( t(i(k),2) + 
     &               p(k)*( 2.*t(i(k),3) + 
     &                      3.*t(i(k),4)*p(k) ) ) / r(k)**2
        pact=pact+t(i(k),1)+p(k)*(t(i(k),2)+p(k)*(t(i(k),3)
     &                                        +p(k)*t(i(k),4)))
      end do
      return
      end
