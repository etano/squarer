      subroutine pbc(xt,x,dx,r,n,ind,rem,csit,cutrt,off,ndim)
      include 'mach.p'
c calculate the distance between xt and x in periodic boundary
c conditions and table entry with remainder
c we use a vectorized logic with unrolled do-loops

      real*8 xt,x,dx,r,rem,csit,cutrt,t,off
      integer n,ind,ndim,nslices,j,l

      dimension ind(n),xt(ndim),x(ndim,n),dx(ndim,n),r(n)
     +,t(3),rem(n)
      include 'cpbc.cm'

c t is needed in statement function fpbct
      do l=1,ndim
         t(l)=sign(ell(l),xt(l))
      end do

      if(ndim.eq.3) then
         do j=1,n
            dx(1,j)=fpbct(xt(1)-x(1,j),1)
            dx(2,j)=fpbct(xt(2)-x(2,j),2)
            dx(3,j)=fpbct(xt(3)-x(3,j),3)
            r(j)=sqrt(dx(1,j)*dx(1,j)
     &                   +dx(2,j)*dx(2,j)     
     &                   +dx(3,j)*dx(3,j) )
            r(j)=min(r(j),cutrt)
            ind(j)=int(r(j)*csit+off)
            ind(j)=max(ind(j),1)
            rem(j)=(r(j)*csit+off)-ind(j)
         end do

      elseif(ndim.eq.2) then
         do j=1,n
            dx(1,j)=fpbct(xt(1)-x(1,j),1)
            dx(2,j)=fpbct(xt(2)-x(2,j),2)
            r(j)=sqrt(dx(1,j)*dx(1,j)+dx(2,j)*dx(2,j))
            r(j)=min(r(j),cutrt)
            ind(j)=int(r(j)*csit+off)
            ind(j)=max(ind(j),1)
            rem(j)=(r(j)*csit+off)-ind(j)
         end do
      endif

      return
      end
