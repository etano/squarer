      function pairdnode(nf,fdenmat,m1,ndim,ifl,r1,r2)
cjbs 6/1/99 created from dnode to work for pair matrices.
      include 'mach.p'
      integer mnp,nf,m1,ndim,ifl,ipr,i,j,k,ize,ip,jp,ivor,l,which
     &   ,i2,i3
      parameter (mnp=256)
      real*8 pairdnode,fdenmat,c,z,r1,r2,dm,dp,ys, yt,d(3)
      real*8 g2,g1,g3,rmin,d2i,gradtemp(ndim,nf)
c finds the distance to the node using linear extrapolation
c nf = rank of matrix, m1= dimension parameter
c ifl=0 use straight Newton, ifl=1 search for P
      logical la(mnp)
      dimension c(mnp,mnp),ize(mnp),ys(mnp),yt(mnp),dm(mnp),dp(mnp)
     +,ivor(mnp)
      dimension fdenmat(m1,ndim+2),ip(mnp),jp(mnp)
     +         ,r1(ndim,nf),r2(ndim,nf)
      save ipr
      data ipr/0/
      if(nf.gt.mnp)then
         write (*,*)' nf>mnp in fnlib/pairdnode ',nf,mnp
         stop
      endif

c find best permutation
      if(ifl.ne.0) then
c optimize over permutation space
      do 50 i=1,nf
      do 50 j=1,nf
        k=j+(i-1)*nf
50      c(j,i)=-log(max(abs(fdenmat(k,1)),1.d-30))
        call lsapr(nf,c,z,ize,ip,dm,dp,ys,yt,ivor,la,mnp)
        !also calculate inverse permutation (do indexing calc here, too)
        do i=1,nf    
          jp(ip(i))=(i-1)*nf
        end do
      else
        !Default to the identity permutation.
        do i=1,nf
          ip(i)=i
          jp(i)=(i-1)*nf
        enddo
      endif

      d(2)=0.d0
      d(3)=0.d0
      !contribution from species1
      do i=1,nf 
        k=(i-1)*nf
        do l=1,ndim 
          g3=0.d0
          do j=1,nf
            g3=g3+ fdenmat(k+j,l+1)*fdenmat(k+j,ndim+2)
          enddo
          ! subtract gradient coming from the diagonal contribution
          ! i.e. take grad of (ln(d) - ln(diag(d)))
          g2=g3-fdenmat(k+ip(i),l+1)/fdenmat(k+ip(i),1)
          d(2)=d(2)+g2**2
          d(3)=d(3)+g3**2
        enddo
      enddo
      !contribution from species2
      do i=1,nf 
        do l=1,ndim 
          g3=0.d0
          do j=1,nf
            k=(j-1)*nf
            g3=g3+ fdenmat(k+i,l+1)*fdenmat(k+i,ndim+2)
          enddo
          ! subtract gradient coming from the diagonal contribution
          ! i.e. take grad of (ln(d) - ln(diag(d)))
          g2=g3-fdenmat(i+jp(i),l+1)/fdenmat(i+jp(i),1)
          d(2)=d(2)+g2**2
          d(3)=d(3)+g3**2
        enddo
      enddo
 
      d(2)=1.d0/sqrt(d(2))
      d(3)=1.d0/sqrt(d(3))
c d(1) is a rigorous upper bound to nodal distance =min(rij)/sqrt(2)
      d(1)=min(rmin(r1,ndim,nf,gradtemp),rmin(r2,ndim,nf,gradtemp))

      i2=2                    !first take minimum between d1 and (d2 and d3)
      if(d(1).lt.d(2)) i2=1
      i3=3
      if(d(1).lt.d(3)) i3=1

      which=i3                !now take maximum of resulting (d2 / d3)
      if(ifl.ne.0.and.d(i2).gt.d(i3)) which=i2

      if(which.lt.3) then !now move data into gradn and dnode
         pairdnode=d(which)
      else
         pairdnode=d(3)
      endif

      return
      end
