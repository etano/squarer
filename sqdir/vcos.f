      subroutine vcos(p,ldim,dc,nc,lu)
      implicit none
      integer l,i,lu,ldim,nc
      real*8 p(0:ldim,nc),dc
      do 1 l=0,lu
       do 1 i=1,nc
1     p(l,i)=cos(dc*(i-1)*l)
      return
      end
