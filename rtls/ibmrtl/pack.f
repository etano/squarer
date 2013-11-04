      subroutine pack(p,nbits,u,nw)
      include 'mach.p'
      integer p(nw),u(nw)
      integer k,ksh,j,nbits,nw,mask
      mask=2**nbits-1
      k=0
      ksh=32
      do 100 j=1,nw
      if(ksh.ge.32) then
        ksh=0
        k=k+1
        p(k)=0
      endif
      p(k)=or( p(k), lshift( and(u(j),mask) ,ksh) )
      ksh=ksh+nbits
100   continue
      return
      end
