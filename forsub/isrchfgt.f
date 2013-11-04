      function isrchfgt(n,v,inc,x)
      include 'mach.p'
c return index of smallest element of ordered vector > x;
c n+1 is returned in case of no success;
      integer n,inc,isrchfgt,ilo,ihi,mid,k
      real*8  v(n),x
      if(inc.ne.1)stop
       if(v(1).gt.x) then
         isrchfgt=1
         return
       endif
      if(v(n).le.x) then
        isrchfgt=n+1
        return
      endif
      ilo=1
      ihi=n
      do 10 k=1,n
        if(ihi.eq.ilo+1) then
          isrchfgt=ihi
          return
        endif
        mid=(ilo+ihi)/2
        if(v(mid).gt.x) then
          ihi=mid
        else
          ilo=mid
        endif
   10 continue
      return
      end
