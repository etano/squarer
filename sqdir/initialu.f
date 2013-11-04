      subroutine initialu(nx,ml,mx2,rv,pott,tau,vint,uold,md)
c this initializes u using pot and tau with Eq. 4.13
      real*8 tau
      integer  nx,ml,mx2,ii,i,l,j,md

      real*8 vint(nx),rv(nx),pott(nx),uold(0:ml,mx2,md)
      write (*,*)' entering initialu '

c vint(i) will be the integral from rv(mx) to rv(i) of v(r)
      vint(nx)=0.d0
      do ii=nx-1,1,-1
        vint(ii)=vint(ii+1)+.5d0*(rv(ii+1)-rv(ii))*(pott(ii+1)+pott(ii))
      enddo
 
c initialize uold; first the diagonal
      do 8 i=1,nx
        do 8 l=0,ml
         uold(l,isym(i,i),1)=tau*pott(i)
         if(md.ge.2)uold(l,isym(i,i),2)=pott(i)
         if(md.ge.3)uold(l,isym(i,i),3)=0.d0
8     continue
 
      do 6 i=1,nx
        do 6 j=1,i-1
c initial density matrix is the integral from ri to rj of potential
         uu=-(vint(i)-vint(j))/(rv(i)-rv(j))
         do 6 l=0,ml
          uold(l,isym(i,j),1)=uu*tau 
          if(md.ge.2)uold(l,isym(i,j),2)=uu  
6         if(md.ge.3)uold(l,isym(i,j),3)=0.d0
      return
      end
