      subroutine unpack(p,nbits,u,nw)
      include 'mach.p'
      integer k,ksh,j,nbits,nw,mask
      integer p(nw),u(nw)
      mask=2**nbits-1
      k=0
      ksh= 32
      do 100 j=1,nw
      if(ksh.ge.32) then
        ksh=0
        k=k+1
c       p(k)=0
      endif
      u(j)=and(mask,rshift(p(k),ksh) )
      ksh=ksh+nbits
100   continue
      return
      end
