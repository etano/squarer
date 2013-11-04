      subroutine readseln(x,mntsc,morder,mterms,mlevels,nterms,nlvld
     +,tau,lptable,ifile,lnew,ifzero)
c read new style tables
c x is where the table goes
c  mntsc,morder,mterms,mlevels are dimensioning parameter (input)
c nterms are the number of terms in approximation (input)
c nlevels = number of levels read ; nlvld=min(nlvld,nlevels)
c tau is the lowest imaginary time needed (input)
c lptable is size of actual table for spline (input)
c ifile = 1,2,3,4,5 for u udot drift and covariance,w (input)
      include 'mach.p'
      real*8 x,tau,dummy,taulow,tauhigh,xk,rlread
      integer mnprm,mntsc,morder,mterms,mlevels,nterms,nlevels,lptable
      integer ifile,idim,k,j,i,n,ipickoff,intread,ksk,l,k0,nlvld,lnew
      dimension x(mntsc,morder,mterms,mlevels),idim(3)
      parameter (mnprm=12)
      character  p(mnprm)*28
      logical ifzero
c first find appropriate table 
2     idim(1)=-mntsc
      idim(2)=0
      idim(3)=0
      call checkrnk(21,3,idim)
      if(idim(1).ne.lptable) then
       write (*,*) ' in readseln : idim(1).ne.lptable '
       stop
      end if
      if(idim(2).ne.nterms) then 
       call findbegn(21)
       read (21,*)(((dummy,i=1,idim(1)),j=1,idim(2)),k=1,idim(3))
       go to 2
      endif

1     if(ipickoff(21,p,n,mnprm).ne.0.or.n.le.0) then
           write (*,*)' data file empty ',n
           stop
      endif
      if(p(1).eq.'BEGIN') then
        if(intread(p(n)).ne.ifile) then
          read (21,*)(((dummy,i=1,idim(1)),j=1,idim(2)),k=1,idim(3))
          go to 2
        else
          go to 10
        endif
      endif
c find tau
      if(p(1).eq.'GRID'.and.intread(p(2)).eq.3) then
        if(p(3).ne.'LOG') then
          write (*,*) ' in readseln : p(3).ne.LOG'
          stop
        end if
        taulow=rlread(p(4))
        tauhigh=rlread(p(5))
        xk=log(tau/taulow)/log(2.)
        ksk=(xk+.001)
        if(abs(xk-ksk).gt..01) then
          write (*,*) ' taugrid on table does not match '
     +       ,tau,taulow,tauhigh,xk
          stop
        endif
      endif

      go to 1

10    continue
      nlevels=0
      k0=1

      do j=1,idim(3)
         do l=1,nterms
            if(k0.gt.mlevels) then
               read (21,*) (dummy,i=1,idim(1))
            else  
               read (21,*) (x(i,1,l,k0) ,i=1,idim(1))
c zero at edge of grid?
               if(ifzero) x(lnew,1,l,k0)=0.d0
               if(j.gt.ksk.and.morder.eq.4)  
     +              call spline(x(1,1,l,k0),lnew,mntsc,1)
            endif
         enddo
         if(j.gt.ksk) nlevels=min(mlevels,nlevels+1)
         if(nlevels.ge.1) k0=k0+1
      enddo

      write (6,*) ' number records skipped ',ksk
     +  ,' levels defined ',nlevels,' terms read ',nterms
      call mcheck(nlevels,mlevels,'nlevels','mlevels','readseln')
      if(nlevels.le.0) then
       write (*,*) ' in readseln : nlevels.le.0'
       stop
      end if
      nlvld=min(nlvld,nlevels)
      return
      end
