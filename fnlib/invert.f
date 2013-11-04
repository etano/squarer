       subroutine invert(n,f,fv,det)
       include 'mach.p'

       integer mnss,n,j,i,ipvt
       real*8 f,fv, det,rcond,work

       parameter (mnss=50)
       dimension f(n,n),fv(n,n),ipvt(mnss),work(mnss),det(2)
       if(n.gt.mnss.or.n.le.0) then
           write(*,*)' mnss too small in invert ',mnss,n
           stop
       endif
c      first transpose
       do j=1,n
       do i=1,n
         fv(j,i)=f(i,j)
       enddo
       enddo
c     call dgefa(fv,n,n,ipvt,info)
c     if(info.gt.0) then
c       write (*,*) ' matrix singular in invert ',info
c       stop
c     endif
 
c estimate the condition number
      call sgeco(fv,n,n,ipvt,rcond,work)
c     write (6,*) ' condition number of dm ',rcond
      if(rcond.lt.1.d-10) then
          write (*,*)' matrix almost singular in invert ',rcond
      endif
      call sgedi(fv,n,n,ipvt,det,work,11)
c sign(det(1))*exp(det(2)) is the value of the determinant
!     det(2)=log(abs(det(1)))+det(2)*log(10.)
      det(2)=log(abs(det(1)))+det(2)*2.302585092994045d0
      return
      end
