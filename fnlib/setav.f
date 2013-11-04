      subroutine setav(ladd,qnamet,aunitt,n)
      include 'mach.p'
      integer ladd,n,iflag,i,k
      include 'caver.cm'
      character  qnamet*(*),aunitt*(*)
      save iflag
       data iflag /0/
       if (iflag.eq.0) then
          naver=0
          iflag=1
       endif
c this subroutine stores qaname aunit and anorm into the arrays in common
c  and increments naver
c go thru previous averages and check for a match
      do 30 i=1,maver-n+1
      do 31 k=1,n
      if(qaname(i+k-1).ne.qnamet)  go to 30
31    continue
c we have found a match with previous averages
      ladd=i
      do 32 k=1,n
      avtemp(i+k-1)=0.0
32    anorm(i+k-1)=0.0
      return
30    continue
3     continue
      ladd=naver+1
      do i=1,n
      naver=naver+1
      if(naver.gt.maver) then
        write (*,*)'  danger too many averages in setav ',naver,maver
        stop
       endif
c zero initial temporary bins
      avtemp(naver)=0.0d0
      anorm(naver)=0.0d0
c store name and units and  zero cumlative bins
      qaname(naver)=qnamet
      aunits(naver)=aunitt
      avsum(naver)=0.0d0
      ablock(naver)=0.0d0
      avsq(naver)=0.d0
      ansum(naver)=0.0d0
      write (63,'(i5,2a18)') 
     .   naver,qaname(naver),aunits(naver)
      enddo
      return
      end
