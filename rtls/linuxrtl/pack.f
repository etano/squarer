      subroutine pack(p,nbits,u,nw)
      include 'mach.p'
      integer nbits,nw,mask,k,ksh,j
      integer p(nw),u(nw)
      mask=2**nbits-1
      k=0
      ksh=32
      do 100 j=1,nw
      if(ksh.ge.32) then
        ksh=0
        k=k+1
        p(k)=0
      endif
!     p(k)= p(k) .or. ( (u(j).and.mask)*2**ksh)
      p(k)=u(j)
      ksh=ksh+nbits
100   continue
      return
      end
