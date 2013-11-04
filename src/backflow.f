      subroutine backflow( uold,diag,rv,lpos,fitdm,chis,norder,
     + ndim, ml,mx2,mx,morder2,morder,nx,rmin,rmax,tau,ehbs2m,bothsz,md)
      implicit none
c transforms uold from partial waves to polynomial method.
      integer mml,ml,mx2,mx,mxx,nx,lpos(nx),lul,nc,ndata,md
     +,morder2,mc,mdata,mo2,ndim,mo,id,korder,k,ii,j,ll,k1
     +,norder,morder,i,i1,ix2,ix3,ix4,isym,l,ktrm,ix,ix1,ic,mdatam
c local parameters: mxx=no of grid points, mml=number of grid points
      parameter (mxx=300,mml=128,mc=30,mdata=10000,mo=10
     +      ,mo2=mo*(mo+3)/2)
       integer nkt(0:mo)
c if bothsz=true we do a 2d fit in s and z. else a 1d fit in s only.
       logical bothsz

      real*8 uold(0:ml,mx2,md),drt(mxx),di(0:mml,3),free(0:mml,2),rv(nx)
     +,f(0:mml,3),uend(3),chis(mx,0:morder,md),ua(3),tau
     +,udata(mdata,mo2),b(mdata,3),v(mo2,mo2),w(mo2)
     +,utemp(mdata,mo2),diag(mx,3)
     +,fitdm(mx,morder2,morder,3),cl(0:mml),u(3)
     +,rvmax,rvmin,rmax,rmin,sc,cutfit,thetam,cmin,dc,dt,u2(3),err(3)
     +,p11,ac(0:mml,mc),ehbs2m,star(3)
     +,t(mo2),ufit(mo2),chimax(0:mo,3),r1,z2,cut,r2,zn,wt,chisq,thresh
     +,swc,wc,ro,rozero,s2,c,pf,wtinv,cll,cns
     +,a1,a2,a3,a4,rhocut,sdot,wmax,sum,rozerod

      if(nx.gt.mxx)stop
      if(norder.gt.mo) stop
      wtinv=.25d0/(tau*ehbs2m)
      pf=(12.5663706143591d0*tau*ehbs2m)**(-.5d0*ndim)
c initialize partial wave components
      call pwc(ndim,mml,cl)
c these are the appended max and min of our grid.
      rvmax=min(2.d0*rv(nx)-rv(nx-2),rmax)
      rvmin=max(2.d0*rv(1)-rv(3),rmin) 
c drt is derivative of the grid
      do i=2,nx-1
          drt(i)=.5d0*abs(rv(i+1)-rv(i-1))
      enddo
      drt(1)=abs(rv(2)-rv(1))
      drt(nx)=abs(rv(nx)-rv(nx-1))
c we cut off fit at 10**-4=e**-7 of peak value
      cutfit=-4*ehbs2m*tau*log(1.d-4)
c integrate to 10^-20 of peak value
        cut=sqrt(-2*ehbs2m*tau*log(1.d-20))
c pointers for fitting expansion
        do  k=0,norder
        if(bothsz) then
           nkt(k)=k*(k+3)/2
        else
           nkt(k)=k
        endif
        enddo
c maximum errors
        do id=1,md
        err(id)=0.
        do k=0,norder
          chimax(k,id)=0.d0
        enddo
        enddo
c scratch space needed
        mdatam=0

c loop over q
      do 600 i=1,nx
c maximum l considered-computed in routine diagonal.
          lul=lpos(i)
          if(lul.gt.mml)stop
c determine maximum value of cos(theta)
          sc=sqrt(cutfit)
          if(sc.ge.rv(i)) then
           thetam=3.14159265d0
           cmin=-1.d0
          else
           thetam=asin(sc/rv(i))
           cmin=cos(thetam)
          endif

c nc is the number of angular integration points
          nc=min(mc,max(1+30*int(1.-cmin),11))
c this  is spacing in cos(theta) (ndim=3) or theta (ndim=2)
          dc=(1.-cmin)/float(nc-1)
          dt=thetam/float(nc-1)
          if(ndim.eq.1)nc=1
c number of fitting points 
          ndata=0
          u2(1)=0.d0
          u2(2)=0.d0
          u2(3)=0.d0
          p11=0.d0

c compute legendre polynomials or cosines
          if(ndim.eq.3)call vlegp(ac,mml,dc,nc,lul)
          if(ndim.eq.2)call vcos(ac,mml,dt,nc,lul)

c now loop over r1
          do 602 i1=1,nx
           r1=rv(i1)
           r2=2.d0*rv(i)-r1
  
c consider (r1,r2) as a fitting point if it is reasonable
           if(abs(r1-r2).le.cut.and.r2.gt.rvmin.and.r2.lt.rvmax) then
           z2=(r1-r2)**2
c calculate powers of z**2
            if(bothsz)then
             zn=1.d0
             do k=1,norder
              zn=zn*z2
              t(nkt(k))=zn
             enddo
            endif
c this is the weight factor 
            wt=(r1*r2)**(ndim-1)*drt(i1)
c now set up density matrix values at r2
            call interp(r2,ix,a1,a2,a3,a4)
            ix1=isym(ix-1,i1)
            ix2=isym(ix,i1)
            ix3=isym(ix+1,i1)
            ix4=isym(ix+2,i1)
c get free particle density matrix 
            call fdenmt(r2,r1,free,tau*ehbs2m,mml,ehbs2m,lul-1,ndim)
c determine interacting part
            do l=0,lul
             di(l,1)=exp(-a1*uold(l,ix1,1)-a2*uold(l,ix2,1)
     +              -a3*uold(l,ix3,1) -a4*uold(l,ix4,1) )
             if(md.ge.2) then
             di(l,2)=a1*uold(l,ix1,2)+a2*uold(l,ix2,2)+a3*uold(l,ix3,2)
     +              +a4*uold(l,ix4,2)
             if(md.ge.3) then
             di(l,3)=a1*uold(l,ix1,3)+a2*uold(l,ix2,3)+a3*uold(l,ix3,3)
     +              +a4*uold(l,ix4,3)
             endif
             endif
            enddo

            do id=1,md
            star(id)=di(lul,id)
            enddo

            cns=(r1*r2)**(-.5d0*(ndim-1))
            do  l=0,lul-1
             cll=cl(l)*cns*free(l,1)
             f(l,1)=cll*(di(l,1)-star(1))
             do id=2,md
               f(l,id)=f(l,1)*free(l,2) 
     +             +cll*(-di(l,id)*di(l,1)+star(id)*star(1))
             enddo
            enddo

c this is the end point value
            do id=1,md
            uend(id)=.5d0*(diag(i1,id)+a1*diag(ix-1,id)+a2*diag(ix,id)
     +               +a3*diag(ix+1,id)+a4*diag(ix+2,id))
            enddo

c we will cut off fitting for rho<rhocut
            rhocut=1.d-4*pf*exp(-uend(1) ) 
c now integrate over cos(theta)
            do 610 ic=1,nc
             c=1.d0-dc*(ic-1)
             if(ndim.eq.2)c=cos(dt*(ic-1))
             s2=r1*r1+r2*r2-2.d0*r1*r2*c

             if(s2.le.cutfit) then

c define independant t variables by recursion
              t(1)=s2
              do k=2,norder
               if(bothsz) then
                do j=1,k
                 t(nkt(k-1)+j)=t(nkt(k-1)+j-k)*s2
                enddo
               else
                t(k)=t(k-1)*s2
               endif
              enddo

              rozero=pf*exp(-wtinv*s2)
              ro=rozero*star(1)+sdot(lul,f(0,1),1,ac(0,ic),1)
             if(ro.gt.rhocut ) then
               u(1)=-log(ro/rozero)
               rozerod=(-.5d0*ndim+wtinv*s2)/tau
               do id=2,md
               u(id)=(rozero*star(1)*(star(id)-rozerod)
     +          -sdot(lul,f(0,id),1,ac(0,ic),1))/ro+rozerod
               enddo
c subtract end-point value
               do id=1,md
                ua(id)=u(id)-uend(id)
               enddo
               if(ic.eq.1.and.i.eq.i1) then
                  do id=1,md
                    err(id)=max(err(id),abs(ua(id)) )
                  enddo
               endif

c add to least squares arrays
               ndata=ndata+1
               if(ndata.gt.mdata) then
                write (*,*)' mdata too small ',mdata
                stop
               endif

               wc=wt*ro
               if(abs(c).gt..999)wc=.5*wc
               swc=sqrt(wc)
               do k1=1,nkt(norder)
                udata(ndata,k1)=swc*t(k1)
               enddo
 
               p11=p11+wc
               do id=1,md
                b(ndata,id)=swc*ua(id)
                u2(id)=u2(id)+wc*ua(id)**2
               enddo

              endif
             endif
610         continue

           endif
602       continue

          mdatam=max0(mdatam,ndata)
c now do least squares fit
          do 700 korder=0,norder
           do 700 id=1,md
 
            if(korder.eq.0) then

             if(u2(id)*p11.le.0.d0) then
              if(id.eq.1)write (*,*) i, ' p11=0 ',p11,rv(i)
              chis(i,0,id)=0.d0
             else
              chis(i,0,id)=sqrt(u2(id)/p11)
             endif

            else
            ktrm=nkt(korder)

             if(id.eq.1) then
c do this once since udata is independant of id
              do 718  k1=1,ktrm
              do 718 ll=1,ndata
718           utemp(ll,k1)=udata(ll,k1)
c now set up rhs and lhs and solve linear equations from NR
              call svdcmp(utemp,ndata,ktrm,mdata,mo2,w,v)
              wmax=0.d0
              do j=1,ktrm
                wmax=max(w(j),wmax)
              enddo
              thresh=1.d-8*wmax
              do j=1,ktrm
               if(w(j).lt.thresh) w(j)=0.d0
              enddo
             endif
           call svbksb(utemp,w,v,ndata,ktrm,mdata,mo2,b(1,id),ufit)
           chisq=0.d0
           do ii=1,ndata
               sum=-b(ii,id)
               do k=1,ktrm
                 sum=sum+ufit(k)*udata(ii,k)
               enddo
               chisq=chisq+sum**2
           enddo
           if(chisq*p11.gt.0.d0)then
               chis(i,korder,id)=sqrt(chisq/p11)
           else
               do k=1,ktrm
                 ufit(k)=0.d0
               enddo
               chis(i,korder,id)=0.d0
           endif
 

c save data. suppress off diagonal elements when diagonal is large
           do k=1,ktrm
              if(diag(i,1).lt.25.) then
                fitdm(i,k,korder,id)=ufit(k)
              else
                fitdm(i,k,korder,id)=0.0
              endif
            enddo
           endif
          chimax(korder,id)=max(chimax(korder,id),exp(-diag(i,1))*
     +    chis(i,korder,id))

700   continue
600   continue
      if(err(1).gt.1.e-4)write (*,*)' diagonal errors ',(err(id),id=1,3)
      write (*,*)' mdatam=',mdatam,'/',mdata
      do korder=0,norder
      write (*,800) korder,(chimax(korder,id),id=1,md)
      enddo
800   format(' order =',i5,' chimax ',3e12.5)

      return
      end
