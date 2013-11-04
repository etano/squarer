      subroutine fillk(rknorm,wtk,kmult,nshlls,nshex,ndim,vol)
      implicit none
c sets up rknorm & wtk for fitp by extending the grid to larger k's
      integer nshlls,kmult(0:nshlls),nshex,ndim,k
      real*8 rknorm(0:nshex),wtk(0:nshex),dk,cutk,pi,con,vdown,vup,vol
!     write (*,*)' fillk', nshlls,nshex,ndim,vol

      wtk(0)=1.d0 ! weight of a given k vector
      do k=1,nshlls
        wtk(k)=2.d0*(kmult(k)-kmult(k-1))
      enddo
c fill in region above cut  with a uniform grid
      dk=rknorm(1)*.5d0
      cutk=rknorm(nshlls)
      pi=3.1415926535d0
      con=vol*(ndim-1)/(ndim*(2.d0*pi)**(ndim-1))
      vdown=1+2*kmult(nshlls)

      do k=nshlls+1,nshex
           vup=con*(cutk+dk*(k-nshlls))**ndim
           rknorm(k)=cutk+dk*(k-nshlls-.5d0)
           wtk(k)=vup-vdown
           vdown=vup
      enddo
      end
