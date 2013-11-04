      subroutine pmult(c,n,np,pn)
      implicit none
      integer k,np,i,n
      real*8 pn,c(n+np)
c multiply by (-r+1)**np
      do k=1,np
           do i=1,n+np-k
              c(i)=c(i+1)-c(i)
           enddo
      enddo
c multiply by pn
      do i=1,n
         c(i)=pn*c(i)
      enddo
      return
      end
