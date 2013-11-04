      subroutine deltar(xt,x,dx,r,n,ndim,nprojn,projn)
      include 'mach.p'
c calculate the distance in pbc
      real*8 xt,x,dx,r,t
      integer i,l,j,nprojn,projn,ndim,nz,n,nslices
      dimension xt(ndim),x(ndim,n),dx(ndim,n),r(n),t(3)
     &  ,projn(2,nprojn)
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
         end do

      elseif(ndim.eq.2) then
         do j=1,n
            dx(1,j)=fpbct(xt(1)-x(1,j),1)
            dx(2,j)=fpbct(xt(2)-x(2,j),2)
            r(j)=sqrt(dx(1,j)*dx(1,j)+dx(2,j)*dx(2,j))
         end do
      endif

c now zero out some components and recompute r
      do j=1,nprojn
      nz=projn(1,j)
      if(nz.le.n) then
      dx(projn(2,j),nz)=0.
c now reform r
      if(ndim.eq.3)r(nz)=sqrt(dx(1,nz)**2+dx(2,nz)**2+dx(3,nz)**2)
      if(ndim.eq.2)r(nz)=sqrt(dx(1,nz)**2+dx(2,nz)**2)
      elseif(nz.eq.n+1)then
        do i=1,n
        dx(projn(2,j),i)=0.
c now reform r
        if(ndim.eq.3)r(i)=sqrt(dx(1,i)**2+dx(2,i)**2+dx(3,i)**2)
        if(ndim.eq.2)r(i)=sqrt(dx(1,i)**2+dx(2,i)**2)
        enddo
      endif
      enddo
c before we do square root should we make a sparse list? of r and dx?

      end

