      subroutine abc(xt,x,n,ind,rem,r,csit,cutrt,off,ndim)
      include 'mach.p'
c calculate the distance between xt and x in periodic boundary
c conditions and table entry with remainder
c we use a vectorized logic with unrolled do-loops
c

      real*8 xt,x,rem,r,csit,cutrt,off,aa
      integer n,ind,ndim,nslices,j

      dimension ind(n),xt(ndim),x(ndim,n),rem(n),r(n)
      include 'cpbc.cm'

      if(ndim.eq.3) then
        do j=1,n
          r(j)=sqrt(fabc(xt(1)-x(1,j),1)
     &             +fabc(xt(2)-x(2,j),2)
     &             +fabc(xt(3)-x(3,j),3))
          r(j)  = min(r(j),cutrt)
          aa    = r(j)*csit+off
          ind(j)= aa
          ind(j)=max(ind(j),1)
          rem(j)=aa-ind(j)
        end do

      elseif(ndim.eq.2) then
        do j=1,n
          r(j)=sqrt(fabc(xt(1)-x(1,j),1)
     &             +fabc(xt(2)-x(2,j),2))
          r(j)  = min(r(j),cutrt)
          aa    = r(j)*csit+off
          ind(j)= aa
          ind(j)=max(ind(j),1)
          rem(j)=aa-ind(j)
        end do
      endif

      return
      end
