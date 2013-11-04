c this is a version from ITP never completed?
      subroutine *goffd(m,i,p,t,n,s,z,norder,mnp,morder,idim,
     +,us,uz,uq)
      include 'mach.p'
c we assume the ordering of terms from squarer
c idim=1 1 d tables , idim=2 full tables.
c us,uz and uq are the returned values of the gradients.

      real*8 p(m),us(mnp),uz(mnp),uq(mnp)
     $, t(n,morder,nterm),s(mnp),z(mnp),var(mnp,0:mno,0:mno)
      integer i(m),n,nterm,mnp,morder,nkt,k,l,j,m,idim

      do j=1,m
       us(j)=0
       uz(j)=0
       uq(j)=0
       var(j,0,0) = 1.
      endif

c now go through the terms in saved order
      nkt=0
      do l=1,norder
      lmin=0
      if(idim.eq.1)ldim=l
      do ls=l,lmin,-1
      lz=l-ls
      nkt=nkt+1

      call tabl(m,i,p,t(1,1,nkt),n,f)
c derivative wrt s
      if(ls.gt.0) then
       do j=1,m
       us(j)=us(j)+ls*f(j)*var(j,ls-1,lz)
       enddo
      endif

c derivative wrt s
      if(lz.gt.0) then
       do j=1,m
       uz(j)=uz(j)+lz*f(j)*var(j,ls,lz-1)
       enddo
      endif

c derivative wrt q
      call gradtabl(m,i,p,t(1,1,nkt),f,n,dxdq,q)

      if(lz.gt.0) then
      do j=1,m
        var(j,ls,lz)=var(j,ls,lz-1)*z(j)
        uq(j)=uq(j)+f(j)*q(j)*var(j,ls,lz)
      enddo

      else
      do j=1,m
        var(j,ls,0)=var(j,ls-1,0)*s(j)
        uq(j)=uq(j)+f(j)*q(j)*var(j,ls,lz)
      enddo
      endif

      enddo
      enddo
      return
      end
