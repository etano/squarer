      subroutine vandvir(m,i,p,t,n,v,vir,rdpdr)
      include 'mach.p'
c given indices i and remainders p finds sum over j of t(rj)
c and r dt/dr 

      real*8 p,t,v,vir,rdpdr
      integer m,i,n,j

      dimension i(m),p(m),t(n,4),rdpdr(m)
      v=0.d0
      vir=0.d0
      do j=1,m
        vir = vir + rdpdr(j)*(t(i(j),2) 
     &                              + p(j)*(2.d0*t(i(j),3) 
     &                                      + 3.d0*p(j)*t(i(j),4)))
         v = v + t(i(j),1) 
     &         + p(j)*(t(i(j),2)
     &                 + p(j)*(t(i(j),3)
     &                         + p(j)*t(i(j),4)))
      end do
      return
      end
