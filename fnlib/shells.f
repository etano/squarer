      subroutine shells(ndim,a,cut,nshlls,rkcomp,rknorm,kmult
     +,nvects,mnkv,mnsh,mdim)
      include 'mach.p'
      integer mnsh
      integer ndim,nshlls,kmult(0:mnsh),nvects,mnkv,npts,l
     +  ,nkspan(3),i,icount(3),nzero,j,kj,jp,mdim
      real*8 a(mdim),rkcomp(mdim,mnkv),rknorm(mnsh)
     +,x(3),cut,c2,rsq,rks,rnow

c computes the vectors x(ndim)=(a(1)*n(1),..,a(ndim)*n(ndim))
c where n(i) are integers and x(1)**2+..+x(ndim)**2.le.cut**2
c the vectors x(i) are stored in rkcomp in
c the order given by the values of their norms.
c  also nshlls gives the number of
c different values of the norms ( the relative square norms
c differ by less than 1.e-5) and cc(lknorm+i) gives
c these nshlls norms and kmult(i) is the last vector
c whose norm is given by knorm(i). hence the total
c number of vectors of magnitude less than cut is kmult(nshlls).
c

      kmult(0)=0
      c2=cut**2
      npts=1
      do 2 l=1,ndim
      nkspan(l)=(0.00001+abs(cut/a(l)))
c range of search is (-nkspan(l),+nkspan(l))
      icount(l)=-nkspan(l)
2     npts=(2*nkspan(l)+1)*npts
      nvects=0
      do 3 i=1,npts
      rsq=0.0d0
c nzero will be the first nonzero entry in icount
      nzero=0
      do 4 l=1,ndim
      if(nzero.eq.0)nzero=icount(l)
      x(l)=icount(l)*a(l)
      rsq=rsq+x(l)**2
      if(rsq.gt.c2) go to 30
4      continue
c we only take half of the vectors and exclude the one at the origin
      if(nzero.le.0) go to 30
      if(nvects.gt.mnkv)
     +call mcheck(nvects,mnkv,'nvects','mnkv','shells')
c we have found a vector
      nvects=nvects+1
c go thru previous vectors. if they have a greater norm move
c them up one slot
      do 6 j=1,nvects-1
      kj=j
       rks=0.d0
      do 66 l=1,ndim
66    rks=rks+rkcomp(l,j)**2
6     if(rks.ge.rsq*1.0001d0) go to 7
      kj=nvects
7     do 8 jp=nvects,kj+1,-1
      do 8 l=1,ndim
8     rkcomp(l,jp)=rkcomp(l,jp-1)
      do 9 l=1,ndim
9     rkcomp(l,kj)=x(l)
c
c
c increase counters with carries for the next vector
30    do 31 l=1,ndim
       icount(l)=icount(l)+1
       if(icount(l).le.nkspan(l)) go to 3
31    icount(l)=-nkspan(l)
3     continue
      nshlls=0
c count number of different norms and find pointers
      rnow=0.d0
      do 51 i=1,nvects
       rsq=0.0d0
       do 514 l=1,ndim
514    rsq=rsq+rkcomp(l,i)**2
      if(rsq-rnow.gt..0001d0*rnow)nshlls=nshlls+1
      if(nshlls.gt.mnsh)
     +         call mcheck(nshlls,mnsh,'nshlls','mnsh','shells')
      rnow=rsq
      rknorm(nshlls)=sqrt(rnow)
51    kmult(nshlls)=i

      end
