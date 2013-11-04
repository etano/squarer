      subroutine writecon(x,ndim,ncomps,ell,fll,mdim)
      include 'mach.p'

      real*8 x,ell
      integer ndim,ncomps,ln,l,i,mdim

      dimension x(mdim,ncomps),ell(mdim)
      character fll*(*)
      ln=index(fll,' ') -1
      open(21,file=fll(1:ln),status='unknown',form='formatted')
      write (21,*) 'RANK ',2,' ',ndim,' ',ncomps
      write (21,22) (ell(l),l=1,ndim)
22    format(' SIZE ',3e18.10)
      write (21,*) 'BEGIN  coordinates'
      do 10 i=1,ncomps
10    write (21,*) (x(l,i),l=1,ndim)
      close(21)
      return
      end
