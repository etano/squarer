      subroutine gvsum(m,i,p,t,var,nterma,nterme,mnp,mntab,morder,mterms
     &,mlevels,iddim,us,uz,uq,rdxdr,r,act,dact)
      include 'mach.p'
      integer mnterms,lmnp,me
      parameter (mnterms=20,lmnp=864,me=2)
cdmc me is the term for the kinetic energy
c we assume the ordering of terms from squarer
c similar to vsum1 but tailored to compute gradients
c idim=1 1-d tables (since s=z) , idim=2 full tables.

      integer m,i(m),nterm,nterma,nterme,mnp,mntab,morder,mterms,mlevels
     &,nkt,k,l,j,iddim,ips(mnterms),ipz(mnterms),ns(mnterms)
     &,nz(mnterms)

      real*8 p(m),t(mntab,morder,mterms,mlevels,3)
     +   ,var(mnp,mterms),utemp(lmnp),rtemp(lmnp),rdxdr(mnp),r(mnp)
     +   ,us(mnp),uq(mnp),uz(mnp),act,dact

         if(m.gt.lmnp)stop
         nterm=max(nterma,nterme)
         if(nterm+1.gt.mterms)stop
c ns,nz and the s and z exponents, ips,ipz are derivative pointers
c   var(nterm+1) will be set to unity
         ns(1)=1
         nz(1)=0
         ips(1)=nterm+1

         ns(2)=0
         nz(2)=1
         ipz(2)=nterm+1

c determine the expansion coefficients
      if(iddim.eq.1) then
c one dimensional case
       ! first copy z**2=var(*,2) down to var(*,1)
        do j=1,m
          var(j,1)=var(j,2)
        enddo
        ! no s dependence in 1d
        ns(1)=0
        nz(1)=1
        ipz(1)=nterm+1
        do l=2,nterm
          ns(l)=0
          nz(l)=l
          ! need exponent z**(2*l-1) in derivative.  This gets
          !   us z**(2l-2).  The other factor of z gets multiplied
          !   back in in "goffd"
          ipz(l)=l-1
          ! build up powers of z**2
          do j=1,m
            var(j,l)=var(j,l-1)*var(j,1)
          enddo
        enddo
      else

c higher dimensional case
       if(nterm.gt.2) then
        nkt=2
        k=1
c go to one higher order. k=ns+nz, nkt= number finished
50      k=k+1
c multiply by s
        do 30 l=1,k
          ns(nkt+l)=ns(nkt+l-k)+1
          ips(nkt+l)=nkt+l-k
          nz(nkt+l)=nz(nkt+l-k)
          ipz(nkt+l)=nkt+l-k-1
        do 30 j=1,m
30      var(j,nkt+l)=var(j,nkt+l-k)*var(j,1)
c multiply by z
          ns(nkt+k+1)=0
          nz(nkt+k+1)=nz(nkt)+1
          ipz(nkt+k+1)=nkt
        do j=1,m
          var(j,nkt+k+1)=var(j,nkt)*var(j,2)
        enddo

        nkt=nkt+k+1
        if(nkt.gt.nterm) then
           write (*,*) 'problem in gvsum ',nterm,nkt
           stop
        elseif(nkt.lt.nterm) then
             go to 50
        endif

        endif
c end of 2-3 dimensional case
      endif
 
      do j=1,m
         rtemp(j)=rdxdr(j)/r(j)
         var(j,nterm+1)=1.d0
         us(j)=0.d0
         uz(j)=0.d0
         uq(j)=0.d0
      enddo

      do l=1,nterm

      if(l.le.nterme) then
        do j=1,m
        dact=dact+(t(i(j),1,l,1,me)
     &                 +p(j)*(t(i(j),2,l,1,me)
     &                        +p(j)*(t(i(j),3,l,1,me)
     &                               +p(j)*t(i(j),4,l,1,me))))*var(j,l)
        enddo
      endif

      if(l.le.nterma) then
        do j=1,m
        uq(j)=uq(j)+rtemp(j)*(t(i(j),2,l,1,1)
     &                        +p(j)*(2*t(i(j),3,l,1,1)
     &                               +3*p(j)*t(i(j),4,l,1,1)))*var(j,l)
        utemp(j)=t(i(j),1,l,1,1)
     &                 +p(j)*(t(i(j),2,l,1,1)
     &                        +p(j)*(t(i(j),3,l,1,1)
     &                               +p(j)*t(i(j),4,l,1,1) ))
       act=act+utemp(j)*var(j,l)

       enddo
c now find s derivatives
       if(ns(l).gt.0) then
        do j=1,m
        us(j)=us(j) + utemp(j)*ns(l)*var(j,ips(l) )
        enddo
       endif
c now find z deriviative
       if(nz(l).gt.0) then
        do j=1,m
        uz(j)=uz(j) + utemp(j)*nz(l)*var(j,ipz(l)  )
        enddo
       endif
      endif
      enddo

      return
      end
