      subroutine echo(lunit,p,n)
      integer n,lunit,nt,i,nc
      character p(n)*(*),line*120
c fill up line with packed p
      if(n.eq.0)return
      nt=1
      do i=1,n
        nc=index(p(i),' ')
        line(nt:nt+nc-1)=p(i)(1:nc)
        nt=nt+nc
      enddo
      if(nt.gt.120)stop
      write (lunit,*) line(1:nt-1)
      return
      end
