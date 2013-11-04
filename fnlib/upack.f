      subroutine upack(x,mnp,ndim,nparts,nslices,el2,ipack,ix,nwp,mcp)
      include 'mach.p'
      real*8 x,el2,c
      integer mnp,ndim,nparts,nslices,ipack,ix,nwp,mcp,l,ind,j,i
      dimension x(ndim,mnp,nslices),ipack(nwp),el2(ndim),ix(nwp),c(3)
!     call unpack(ipack,16,ix,mcp*nwp) 
c now map ix into a real number    
!     do 4 l=1,ndim 
!     c(l)=(2*el2(l))/65536.d0
!     enddo
!     ind=0
!     do  j=1,nslices    
!     do  i=1,nparts 
!     do  l=1,ndim 
!     ind=ind+1   
!     x(l,i,j)=c(l)*(ipack(ind)+.5)-el2(l) 
!     x(l,i,j)=c(l)*(ix(ind)+.5)-el2(l) 
!     enddo
!     enddo
!     enddo
      end 
