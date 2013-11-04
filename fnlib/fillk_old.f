      subroutine fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol
     +         ,mnsh)
      implicit none
      integer mnsh,kmult(0:mnsh),nshlls,nshex1,ndim,k,nshex
c sets up rknorm wtk for fitp routines by extending the grid to larger k's
      real*8 rknorm(mnsh),wtk(mnsh),cutk,pi,con,vdown,vup,dk,sang,argek
     .,vol
      rknorm(1)=0.d0
c this is the weight of a given k vector
      wtk(1)=1.d0
      do k=1,nshlls
        wtk(k+1)=2.d0*(kmult(k)-kmult(k-1))
      enddo
c fill in region between cut and argek with a uniform grid
      dk=rknorm(2)*.5d0
      cutk=rknorm(nshlls)
      nshex=nshlls+(argek-cutk)/dk
      call mcheck(nshex+1,mnsh,'nshex','mnsh','fillk')
      pi=3.1415 92653 58979d0
      con=vol*sang(ndim)/(ndim*(2*pi)**ndim)
      vdown=1+2*kmult(nshlls)
      do k=nshlls+1,nshex
           vup=con*(cutk+dk*(k-nshlls))**ndim
           rknorm(k+1)=cutk+dk*(k-nshlls-.5d0)
           wtk(k+1)=vup-vdown
           vdown=vup
      enddo
      nshex1=nshex+1
      return
      end
