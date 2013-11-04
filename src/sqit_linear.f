      subroutine sqit_linear(nx,nl,ml,mx2,ndim,tau,ehbs2m,rv,rmin,rmax
     +,vint,uold,unew,di,dj,md,tauk,utau)

      implicit none
      integer nx,ml,nl,i,j,k,l,isym,id,indl,mx2,ndim,is,kp,ind,nk
     +,jl,lul,nint,icase,md
      parameter (nk=8)
      real*8 xint(nk),weight(nk),width,cut,etau,rmin,rmax,rv(nx),rr2
     +,uold(0:ml,mx2,md),unew(0:ml,mx2,md),delr,cutl,rbar
     +,rk,wt,ehbs2m,di(0:ml,3),dj(0:ml,3),temp,rinv,tau,up(3),vint(nx)
     +,ektau,tauk,utau(0:ml,mx2,md)
      save xint,weight

c these hermite coefs are from a-s table 25.10
      data xint/.27348104613815d0,.82295144914466d0,1.38025853919888d0
     +,1.95178799091625d0, 2.54620215784748d0,3.17699916197996d0
     +,3.86944 79048 6012d0,4.68873 89393 0582d0/
      data weight/.54737 52050 378d0,0.55244 19573 675d0,.56321 78290 882d0
     +,.58124 72754 009d0,.60973 69582 560d0,.65575 56728 761d0
     +,.73824 56222 777d0, .93687 44928 841d0/

      etau=tau*ehbs2m
      ektau=tauk*ehbs2m
      if(nl.gt.ml) then
        write (*,*) ' di too small in sqit '
        stop
      endif
      width=sqrt(ektau+etau)
c integrate to 10^-20 of peak value
      cut=width*sqrt(-log(1.d-20))
      cutl=-log(1.d-6)

c loop over all i and j values
      do 110 i=1,nx
      do 110 j=i,1,-1
      ind=isym(i,j)
      if(abs(rv(j)-rv(i)).lt.cut) then
           jl=j
           lul=amin0(nl,2+int(sqrt(cutl*rv(i)*rv(j))/width ) )  

           rbar=.5d0*(rv(i)+rv(j))
           do l=0,lul
            unew(l,ind,1)=0.d0
           enddo

           icase=0
           if(rbar-width*xint(nk).le.rmin) icase=1
           if(rbar+width*xint(nk).ge.rmax) icase=icase+2


           if(icase.eq.0) then
c gaussian integration
            nint=nk+nk
           elseif(icase.eq.1) then
c wall is only on the left
            rk=rmin
            delr=0.1d0*width
            nint=int((rbar-rmin+width*xint(nk))/delr)
           elseif(icase.eq.2) then
c wall is only on the right
            rk=rmax
            delr=-0.1d0*width
            nint=int(-(rmax-rbar+width*xint(nk))/delr)
           elseif(icase.eq.3) then
c wall is both on left and right
            rk=rmin
            nint=41
            delr=(rmax-rmin)/(nint+1.d0)
           endif

c now integrate over r
      do k=1,nint

            if(icase.gt.0) then
             rk=rk+delr
             wt=.66666 66666 6666d0*abs(delr)
             if(mod(k,2).eq.1)wt=wt+ wt
            else
             kp=(k-1)/2+1
             is=mod(k,2)
             rk=rbar+width*(2*is-1)*xint(kp)
             wt=weight(kp)*width
            endif
c calculate the density matrix between rk and ri and rk and rj
            call denmat(rk,i,rv(i),uold,ml,nx,di,ektau,ehbs2m,mx2
     +              ,lul,ndim,md)
            call denmat(rk,j,rv(j),utau,ml,nx,dj,etau,ehbs2m,mx2
     +              ,lul,ndim,md)
            do l=0,lul
              unew(l,ind,1)=unew(l,ind,1)+wt*di(l,1)*dj(l,1)
            enddo

      enddo
c now find the free-particle action and divide it out.
           call fdenmt(rv(i),rv(j),di,ektau+etau,ml,ehbs2m,lul,ndim)
           do l=0,lul
            temp=-log(unew(l,ind,1)/di(l,1))
            unew(l,ind,1)=temp
           enddo

c linearly extrapolate in l
           do l=lul+1,nl
            unew(l,ind,1)=2.d0*unew(l-1,ind,1)-unew(l-2,ind,1)
           enddo

      else
c go outwards from the diagonal using classical value instead
           rinv=1.d0/(rv(j)-rv(i))
           up(1)=(tau+tauk)*up(2)
           rr2=rinv*(rv(jl)-rv(i))
           indl=isym(i,jl)
           do l=0,nl
            do id=1,md
            unew(l,ind,id)=up(id)+rr2*unew(l,indl,id)
            enddo
           enddo
      endif

110   continue

c copy new values of density matrix into old array
         do 199 i=1,nx
          do 199 j=1,i
           ind=isym(i,j)
           do 199 l=0,nl
199   uold(l,ind,1)=unew(l,ind,1)
      return
      end
