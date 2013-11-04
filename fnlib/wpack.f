      subroutine wpack(x,nparts,ndim,el2,ix,ipack,nslices,mnp,join,iref
     +,irev,iperm,mcp)
      include 'mach.p'
c call this routine to write packed configurations on file 90

      real*8 x,el2,ci
      integer nparts,ndim,ix,ipack,nslices,mnp,join,iref,irev,iperm
      integer mcp,l,ind,j,i,nwp,nwds

      dimension el2(ndim),x(ndim,mnp,nslices),ix(nparts),ipack(nparts)
     &        , ci(4),iperm(mnp)
c x contains the coordinates in range -el2 to +el2
c nparts and ndim are as in ipwrt
c ix and ipack are work space of size nparts*ndim*nslices and ix/4
c first map into integers in the appropriate range
!     do  l=1,ndim
!       ci(l)=65536/(2*el2(l))
!     enddo
!     ind=0
!     do 2 j=1,nslices
!      do 2 i=1,nparts
!      do 2 l=1,ndim
!      ind=ind+1
!2     ix(ind)=ci(l)*(x(l,i,j)+el2(l))
!     nwp=1+(ind-1)/mcp
!      nwp=ind
!      nwds=mcp*nwp
c now pack ix into ipack with cray library routine
!     call pack(ipack,32,ix,nwds)
c write it to file 90
!     write (90) join,iref,irev,(iperm(i),i=1,nparts)
!    .,(((x(l,i,j),l=1,ndim),i=1,nparts),j=1,nslices)
!     write (90) join,iref,irev,(iperm(i),i=1,nparts),(ipack(i),i=1,nwp)
      do i=1,nparts !straight ascii write statement:
!     do j=1,nslices
      do j=1,1 ! only first slice`
      write (90,'(3f10.4)') (x(l,i,j),l=1,ndim)
      enddo
      enddo
      end
