      subroutine multiup(ratinvo,nord,label,lstart,nc,columns
     +,work,lc,ifrorc,det,ratinvn)
      include 'mach.p'
      real*8 ratinvo,columns, work,det,ratinvn,di
      integer nord,label,lstart,nc,lc,ifrorc,i,j,jp
      dimension ratinvo(nord,nord),label(nc),columns(lc,nord)
     +,work(nord,2),ratinvn(nord,nord)
c updates several rows/columns and puts them into ratinvn

      do 20 i=1,nord
      do 20 j=1,nord
20    ratinvn(i,j)=ratinvo(i,j)
      det=1.
      do 10 i=1,nc
      jp=label(i)-lstart+1
      if(ifrorc.eq.1) then
c column update
        call chgclm(ratinvn,nord,nord,columns(1,i),jp,work,work(1,2),di)
      else
c row update
        call chgrow(ratinvn,nord,nord,columns(1,i),jp,work,work(1,2),di)
      endif
      det=det*di
10    continue
      return
      end
