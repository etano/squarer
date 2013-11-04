      subroutine denmat(r,i,ri,u,ldim,nx,d,et,ehbs2m,mx2,lu,ndim,md)
      implicit none
      integer i,ldim,nx,mx2,lu,ndim,ix,i1,i2,i3,l,isym,i4,md,id
      real*8 r,ri,et,ehbs2m,a1,a2,a3,d2,a4
      real*8 u(0:ldim,mx2,md),d(0:ldim,3)
c finds density matrix at point (r,r(i),et)
c first get free particle density matrix and beta derviative of log(rho_0)
      call fdenmt(r,ri,d,et,ldim,ehbs2m,lu,ndim)
c get quadratic interpolation coefficients
      call interp(r,ix,a1,a2,a3,a4)
      i1=isym(ix-1,i)
      i2=isym(ix,i)
      i3=isym(ix+1,i)
      i4=isym(ix+2,i)
c loop over partial waves, multiplying by free-particle and interpolating
      do l=0,lu
c this is the density matrix
       d(l,1)=d(l,1)*exp(-a1*u(l,i1,1)-a2*u(l,i2,1)-a3*u(l,i3,1) 
     +  -a4*u(l,i4,1))
c this is the beta derivative of the density matrix
       d2=d(l,2)
       do id=2,md
       d(l,id)=d(l,1)*(d2-a1*u(l,i1,id)-a2*u(l,i2,id)-a3*u(l,i3,id)
     +  -a4*u(l,i4,id))
       enddo
      enddo
      return
      end
