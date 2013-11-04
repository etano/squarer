      subroutine readcon(x,ndim,ncomps,ell,filen)
      include 'mach.p'
      real*8 x,ell
      integer ndim,ncomps,ln,idim,i,j
      character filen*(*)
      dimension x(ndim,ncomps),ell(ndim),idim(2)
c reads a configuration file
      ln=index(filen,' ')-1
      open(21,file=filen(1:ln),status='old',form='formatted')
      write (6,*)' beginning read from file ',filen
      idim(1)=ndim
      idim(2)=ncomps
      call checkrnk(21,2,idim)
      call findbegn(21)
      read (21,*) ((x(i,j),i=1,ndim),j=1,ncomps)
      close(21)
      return
      end
