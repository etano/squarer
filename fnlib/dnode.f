      function dnode(nf,fdenmat,m1,ndim,ifl,r,gradn)
cdmc 11/4/96 changed to return unit vector away from node. gradn
      include 'mach.p'
      integer mnp,nf,m1,ndim,ifl,i,j,k,ize,ip,ivor,l,which
     &   ,i2,i3
      parameter (mnp=512)
      real*8 dnode,fdenmat,c,z,r,dm,dp,ys, yt,d(3)
      real*8 g2,g1,g3,rmin,d2i,gradn(ndim,nf),gradtemp(ndim,mnp,2)
c finds the distance to the node using linear extrapolation
c nf = rank of matrix, m1= dimension parameter
c ifl=0 use straight Newton, ifl=1 search for P
      logical la(mnp)
      dimension c(mnp,mnp),ize(mnp),ys(mnp),yt(mnp),dm(mnp),dp(mnp)
     +,ivor(mnp)
      dimension fdenmat(m1,ndim+2),ip(mnp),r(ndim,nf)
      if(nf.gt.mnp) then
         write (*,*)' nf>mnp  in fnlib/dnode ',nf,mnp
         stop
      endif

c find best permutation
      if(ifl.ne.0) then
      d(2)=0.d0
c optimize over permutation space
      do i=1,nf
      do j=1,nf
        k=j+(i-1)*nf
        c(j,i)=-log(max(abs(fdenmat(k,1)),1.d-30))
      enddo
      enddo
      ! solve linear sum assignment problem for best permutation(ip). c=cost
      call lsapr(nf,c,z,ize,ip,dm,dp,ys,yt,ivor,la,mnp)
      endif

 
c     if(ifl.ne.0) then
c     tr=1.
c     do i=1,nf
c        io=(i-1)*nf
c        tr=tr*fdenmat(ip(i)+io,1)
c     enddo
c     if(abs(log(tr)+z).gt..0001) then
c         write (*,*)' problem in laspr ',log(tr),-z,(ip(i),i=1,nf)
c         stop
c     endif
c     endif

      d(3)=0.d0
      do i=1,nf 
      k=(i-1)*nf
      do l=1,ndim 
         g3=0.d0
         do j=1,nf
            g3=g3+ fdenmat(k+j,l+1)*fdenmat(k+j,ndim+2)
         enddo
         gradn(l,i)=g3  
         d(3)=d(3)+g3**2
         if(ifl.ne.0) then
           g2=g3-fdenmat(k+ip(i),l+1)/fdenmat(k+ip(i),1)
!        subtract gradient coming from the diagonal contribution
!        i.e. take grad of (ln(d) - ln(diag(d)))
           gradtemp(l,i,2)=g2
           d(2)=d(2)+g2**2
        endif
      enddo
      enddo
 
c d(1) is a rigorous upper bound to nodal distance =min(rij)/sqrt(2)
      d(1)=rmin(r,ndim,nf,gradtemp)

      d(3)=1.d0/sqrt(d(3))
      i3=3
      if(d(1).lt.d(3)) i3=1
      which=i3

      if(ifl.ne.0) then
          d(2)=1.d0/sqrt(d(2))
          i2=2                    
          if(d(1).lt.d(2)) i2=1
          if(d(i2).gt.d(i3)) which=i2
      endif

      if(which.lt.3) then !now move data into gradn and dnode
         dnode=d(which)
         do l=1,ndim
         do j=1,nf
          gradn(l,j)=gradtemp(l,j,which)
         enddo
         enddo
      else
         dnode=d(3)
      endif
      if(dnode.le.0.0) then
          write (*,*)' problem in dnode ',dnode
          stop
      endif
!     write (*,*) dnode,ifl,which,d(1),d(2),d(3)
!     stop

      return
      end
