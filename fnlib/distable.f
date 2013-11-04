         subroutine distable(r1,r2,n,h,mh,ndim,wi,iflag,cump)
c sets up a matrix h of distances between coordiates r1 and r2
         include 'mach.p'
         integer n,mh,ndim,i,j,iflag,nslices
         real*8 wi,r1(ndim,n),r2(ndim,n),h(mh,n),cump(mh,n),w(3)
         include 'cpbc.cm'

         if(ndim.eq.3) then
          w(1)=wi
          w(2)=wi
          w(3)=wi
          do j=1,n
           do i=1,n
            h(i,j) = exp(-w(1)* fabc(r1(1,i)-r2(1,j),1)
     &                   -w(2)* fabc(r1(2,i)-r2(2,j),2)
     &                   -w(3)* fabc(r1(3,i)-r2(3,j),3) )
          enddo
          enddo
         elseif(ndim.eq.2) then
          w(1)=wi
          w(2)=wi
          do  j=1,n
           do  i=1,n
            h(i,j) = exp(-w(1)* fabc(r1(1,i)-r2(1,j),1)
     &                   -w(2)* fabc(r1(2,i)-r2(2,j),2) )
          enddo
          enddo
         else
           stop
         endif
c if ipf=2 find cumlative distribution
         if(iflag.eq.2) then
           do j=1,n
             cump(1,j)=h(1,j)
           enddo

           do i=2,n
             do j=1,n
               cump(i,j)=cump(i-1,j)+h(i,j)
             enddo
           enddo
         endif

         end
