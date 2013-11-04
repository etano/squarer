      function vsum1(m,i,p,t,n,var,nterm,mnp,morder,idim)
      include 'mach.p'
c we assume the ordering of terms from squarer
c idim=1 1 d tables , idim=2 full tables.
      integer m,i(m),n,nterm,mnp,morder,nkt,k,l,j,idim
      real*8 vsum1,p(m),var(mnp,nterm),t(n,morder,nterm)

c determine the expansion coefficients
      if(idim.eq.1) then

         do l=2,nterm
         do j=1,m
             var(j,l)=var(j,l-1)*var(j,1)
         enddo
         enddo

      else
       if(nterm.gt.2) then
        nkt=2
        k=1

50      k=k+1 ! go to one higher order
        do l=1,k
        do j=1,m
            var(j,nkt+l)=var(j,nkt+l-k)*var(j,1) ! multiply by s**2
         enddo
         enddo
        do j=1,m
          var(j,nkt+k+1)=var(j,nkt)*var(j,2) ! multiply by z**2
        enddo

        nkt=nkt+k+1
        if(nkt.gt.nterm) then
           write (*,*) 'problem in vsum1 ',nterm,nkt
           stop
        elseif(nkt.lt.nterm) then
             go to 50
        endif

        endif
      endif
 
c now do table look-up. for nterm large the loops should be inverted
      vsum1=0.d0
      do l=1,nterm
        do j=1,m
           vsum1=vsum1+(t(i(j),1,l)
     &           +p(j)*(t(i(j),2,l)
     &           +p(j)*(t(i(j),3,l)
     &            +p(j)*t(i(j),4,l)  )))*var(j,l)
      enddo
      enddo
      end
