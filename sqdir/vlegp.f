      subroutine vlegp(p,ldim,dc,nc,lu)
      implicit none
c vector legendre polynomial table
      integer ldim,nc,lu,i,l
      real*8 p(0:ldim,nc),dc,x,a1,a2
      do  i=2,nc
           p(0,i)=1.d0
           p(1,i)=1.d0-dc*(i-1)
      enddo
c use recurrence relations in the increasing l direction
      do l=1,lu-1
       x=1.d0/(l+1)
       a1=x*(l+l+1)
       a2=-x*l
       do i=2,nc
          p(l+1,i)=p(1,i)*a1*p(l,i)+a2*p(l-1,i)
       enddo
      enddo

      do l=0,lu
         p(l,1)=1.d0
      enddo
      return
      end
