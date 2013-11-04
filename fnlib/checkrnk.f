      subroutine checkrnk(iu,nr,idim)
      include 'mach.p'
      integer mnprm,iu,nr,idim, n,i,nread,ipickoff,intread
      parameter (mnprm=12)
      character  p(mnprm)*28
      dimension idim(nr)
c checks a standard data file for rank, unit=iu is assumed to be open
c idim are the number of elements expected if positive
c if idim=0, idim is read in and not checked
c if idim<0 value read in is checked to be .le. than -idim
c stop if eof or empty
1     if(ipickoff(iu,p,n,mnprm).ne.0.or.n.le.0) then
        write (*,*)' data file empty ',n
        stop
      endif

      if(p(1).eq.'RANK') then
c stop if desired rank is greater than read in rank
        if(nr.eq.0) then
          nr=intread(p(2))
        elseif(intread(p(2)).lt.nr) then
          write(*,*)' nrank wrong in checkrnk ',intread(p(2)),nr
          stop
        endif
c now check dimension statements
        do 10 i=1,nr
          if(i.eq.nr.and.nread.eq.0)then
            write(*,*)' last dim read will be finished on eof'
          else
            nread=intread(p(2+i))
            if(nread.lt.0) then
              write (*,*)' array bounds nonsense on ',iu,i,nread
              stop
            endif
            if(idim(i).lt.0.and.nread.gt.-idim(i)) then
              write (*,*)i,' array too small checkrnk ',-idim(i),nread
              stop
            endif
            if(idim(i).le.0) then
              idim(i)=nread
            else 
              if(idim(i).ne.nread) then
                write(*,*)i,' ndim in checkrnk ',nread,idim(i)
                stop
              endif
            endif
          endif
10      continue
        return
      endif
      go to 1
      end
