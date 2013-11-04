      subroutine newup(f,g,c,dn,n,movers,nfirst,nmovers,ndim,ipopt
     +    ,m1,mc,sdet,idebug,ifrorc,r,gradn)
c generates new matrices etc. caused by a change in some columns or rows
c n = rank of matrix; m1 mc ndim are dimensioning parameters
c nmovers = number of changed columns/rows (ifrorc =1 for columns)
c added gradn to find gradient of nodal distance.
      include 'mach.p'

      integer mnss,n,movers,nfirst,nmovers,ndim,ipopt,m1,mc,idebug
      integer ifrorc,l,i,k,ipvt,kk,ic,j,inext,ibv,ll,mmovers
      real*8 f,g,c,dn,sdet,r,dmax,dd,ds,rowi,dot,dt,det,dnode
     +,   gradn(ndim,n)

      parameter (mnss=66,mmovers=5)
      dimension f(m1,ndim+2),g(m1,ndim+2),c(mc,ndim+1),r(ndim,n)
     +,movers(nmovers),det(2),dot(mnss),rowi(mnss),ipvt(mmovers)
     +,ds(mmovers)

      if(n.gt.mnss.or.n.le.0) then
         write (*,*)' mnss too small in newup ',n,mnss
         stop
      endif
      if(nmovers.gt.mmovers)stop

c see if old update is good enough
      if(idebug.ge.2)call checkin(n,f,f(1,ndim+2),66)
      sdet=1.d0
c first copy f into g  (does f=g?)
      do 12 l=1,ndim+2
      do 12 i=1,n*n
12    g(i,l)=f(i,l)

      do 76 k=1,nmovers
76    ipvt(k)=0
c now put in new columns or rows
      do 30 k=1,nmovers

c decide on which column to update next
        dmax=-1.d0
        do 70 kk=1,nmovers
        if(ipvt(kk).eq.0) then
          ic=movers(kk)-nfirst+1

c    find determinant for this update
          dd=0.d0
          if(ifrorc.eq.1) then
          do j=1,n !   inverse*kth column
            dd=dd+g(j+(ic-1)*n,ndim+2)*c(j+(kk-1)*n,1)
          enddo
          else
c this is a bad stride if ifrorc.ne.1 we need gttranspose
          do j=1,n
            dd=dd+g(ic+(j-1)*n,ndim+2)*c(j+(kk-1)*n,1)
          enddo
          endif

          if(abs(dd).gt.dmax) then
                inext=kk
                dmax=abs(dd)
                ds(k)=dd
          endif
        endif
70      continue

        if(abs(ds(k)).le.1.d-7) then
              ibv=1
              go to 50
!             write (*,*)' pivot 0 in newup ',ifrorc,nmovers,k
!             sdet=0.d0
!             return
        endif

      ipvt(inext)=k
      ic=movers(inext)-nfirst+1
c this is the Sherman-Morrison update formula for the inverse
      if(ifrorc.eq.1) then
      do 23 l=1,ndim+1
      do 23 j=1,n
23    g(j+(ic-1)*n,l)=c(j+(inext-1)*n,l)
      call chgclm(g(1,ndim+2),n,c(1+(inext-1)*n,1),ic,rowi,dot,ds(k))
      else
c another bad stride
      do 24 l=1,ndim+1
      do 24 j=1,n
24    g(ic+(j-1)*n,l)=c(j+(inext-1)*n,l)
      call chgrow(g(1,ndim+2),n,c(1+(inext-1)*n,1),ic,rowi,dot,ds(k))
      endif
      sdet=ds(k)*sdet
30    continue

c  check new inverse
      if(idebug.ge.2)call checkin(n,g,g(1,ndim+2),67)
      ibv=0
      do j=1,n
	 dt=-1.d0
	 ll=1+(j-1)*n
	 do i=1,n
	    dt=dt+g(ll,1)*g(ll,ndim+2)
	    ll=ll+1
         end do
         if(abs(dt).gt.1.e-3) then
            write (6,*)' update bad ',dt,j,ifrorc
            write (6,*) ' ds movers ',(ds(k),k=1,nmovers)
     +                   ,(movers(k),k=1,nmovers)
            ibv=1
c           stop
            go to 50
         endif
      enddo
c if flag is set redo inverse
50    continue
      if(ibv.eq.1) then
      call invert(n,g,g(1,ndim+2),det)
      endif
     
      dn=dnode(n,g,m1,ndim,ipopt,r,gradn)
      return
      end
