      subroutine newdet(v,n,label,lstart,nc,columns,dnew,ifrorc)
      include 'mach.p'
c finds the ratio of the new determinant to the old determinant when 
c nc columns are changed.  v is the old inverse matrix (transposed).
c ifrorc =1 means column update otherwise row update
c the columns(rows) to be changed are labeled label(i)-lstart+1
c result goes into dnew

      integer mnmovers,n,label,lstart,nc,ifrorc,i,j,jp,k
      real*8 v,columns,dnew,d

      parameter (mnmovers=3)
      dimension v(n,n),columns(n,nc),d(mnmovers,mnmovers)
     +,label(nc) 

      if(nc.gt.mnmovers.or.nc.le.0) then
        write (*,*) ' nc out of range in newdet ',nc,mnmovers
        stop
      endif

      do 10 i=1,nc
      do 10 j=1,nc

      d(i,j)=0.
      jp=label(j)-lstart+1

      if(jp.le.0.or.jp.gt.n) then
         write(*,*)' labels out of range in newdet ',nc,lstart
     +   ,(label(k),k=1,nc)
         stop
      endif

      if(ifrorc.eq.1) then
        do 20 k=1,n
20      d(i,j)=d(i,j)+columns(k,i)*v(k,jp)
      else
        do 30 k=1,n
30      d(i,j)=d(i,j)+columns(k,i)*v(jp,k)
      endif

10    continue
c     now find determinant of d
      if(nc.eq.1) then
        dnew=d(1,1)
      elseif(nc.eq.2) then
        dnew=d(1,1)*d(2,2)-d(1,2)*d(2,1)
      elseif(nc.eq.3) then
        dnew= d(1,1)*(d(2,2)*d(3,3) - d(3,2)*d(2,3))
     +      + d(1,2)*(d(2,3)*d(3,1) - d(3,3)*d(2,1))
     +      + d(1,3)*(d(2,1)*d(3,2) - d(3,1)*d(2,2))
       endif

       return
       end
