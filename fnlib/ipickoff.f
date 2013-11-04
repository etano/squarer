      function ipickoff(iu,p,n,mn)
      include 'mach.p'
c on entry: iu is unit number for read, mn is the size of array p
c on exit pickoff.eq.1 means eof was encountered
c p contains as characters the parameters, there are n of them
c  blanks or commas separate characters
c line length is maximally 320 and characters after a : or ; are ignored.

      integer line,ipickoff,iu,n, mn,ifact,i,k,ld,ist,ln,ls

      parameter (line=360)
      character data*361, p(mn)*(*)
100   continue
      ifact=0
      n=0 
      i=0
      read (iu,6,END=2) data
      ipickoff=0
6     format(a360)
c starting from end find how long the line is
      do k=line,1,-1
        if(data(k:k).ne.' ') go to 401
      enddo
401   ld=k

      if(iu.le.1)write (6,7)iu,data(1:ld)
7     format(i3,' CMD:',a)
c now strip off comment part
      ls=0
      do k=1,ld
      if(data(k:k).eq.':'.or.data(k:k).eq.';'.or.data(k:k).eq.'!') 
     &  go to 402
       ls=k
      enddo
402   ls=ls+1
      data(ls:ls)=' '

1     i=i+1
      if(i.gt.ls)then
        if(n.gt.0)then
c          write (6,*) ls,ld,n,'=n ',(p(i),i=1,n)
           return
         endif
c read a new record if a blank line is encountered
         go to 100
       endif

      if(data(i:i).eq.' '.or.data(i:i).eq.',')then

       if(ifact.ne.0) then
c finish the word
        ln=i-ist
        if(ln.ge.len(p(n))) then
          write (*,*) 'space for parameters',len(p(n))
     &               ,' not large enough in pickoff',ln
          stop
        endif

        p(n)(1:ln)=data(ist:i-1)
c terminate with  blanks
        do 30 k=ln+1,len(p(n))
30      p(n)(k:k)=' '
        ifact=0
       endif

      else

       if(ifact.eq.0) then
c found a new parameter
         ifact=1
         ist=i
         n=n+1
         if(n.gt.mn) then
           write (*,*)' too many parameters in pickoff',n,mn
           stop
         endif
       endif

      endif
      go to 1

2     ipickoff=1
      return
      end
