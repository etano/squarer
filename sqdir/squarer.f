      program squarer
      implicit none

c nx(mx) are the number of radial grid points
c nl(ml) are the number of partial waves
c norder(morder) is the order of the polynomial fit 
c ntemp(mtemp) are the number of temperatures to be output
      integer mx,nx,ml,nl,norder,morder,ntemp,mtemp,mx2,morder2,mn,md

      parameter (mx=999, ml=99, morder=3, mtemp=12, md=2)
c md= number of derivatives.1: u only, 2:also ub, 3:also uk.
c 
      parameter (mx2=mx*(mx+1)/2, morder2= (morder+3)*morder/2)
c mn is the maximum number of parameters in an input line
      parameter (mn=20)

      real*8 ehbs2m,tlow,taulow,tauhigh,rmin,rmax,hbs2m,rlread,tau
     +,uold(0:ml,mx2,md),unew(0:ml,mx2,md),pott(mx),sum(0:md),wt
     +,vint(mx),rv(mx),fitdm(mx,morder2,morder,md,mtemp),sa,sang
     +,chis(mx,0:morder,md,mtemp)
     +,di(0:ml,3),dj(0:ml,3),diag(mx,md,mtemp),rho,bf(mx,mtemp)

      integer nkt(0:morder),lpos(mx),i,ln,ntypes,ipoint,ires,ipickoff
     +,n,ngrid,intread,k,nsquare,isq,id,ns,j,l,ndim,ifdump,lmax,itemp
     +,ind,isym,ndump,idump(5),isd

      logical ifex,bothsz

      character fname*30,p(mn)*28,grid(mn)*28
      character *100 BUFFER
 
      common /ipp/ipoint(1000)

c set up pointers to store 1/2 of matrix
      ipoint(1)=0
      do i=2,mx
        ipoint(i)=ipoint(i-1)+i-1
      enddo
c in case no potentials are read in
      do i=1,mx
        pott(i)=0.d0
      enddo
      bothsz=.true.

      write (*,*)' Squarer: computes the thermal density matrix'
      write (*,*)' Author: D. Ceperley, University of Illinois. Version
     + of March 10, 1996'
c which file to read??

      call getarg(1,BUFFER)
      read(BUFFER,*) fname

41    write (*,*)' input prefix of files '
c      read (*,99) fname
99    format(a30)
      ln=index(fname,' ')-1
c does this file exist?
      inquire(file=fname(1:ln)//'.dm',exist=ifex)
      if(ifex) go to 40
          write (*,*)' file=',fname(1:ln),'.dm does not exist '
          go to 41
40    open (1,file=fname(1:ln)//'.dm',status='unknown')
c*** read input file dm

c we zero these to ensure proper data is later read in.
      ntypes=0
      ehbs2m=0.d0
      tlow=0.d0
      nx=0
      ifdump=0
c  density matrix is zeroed at rmin and rmax
      rmin=-1.d90
      rmax=1.d90

c read a new input line into p  or stop (n is the number of parameters)
1     ires=ipickoff(1,p,n,mn)
c this jump for zero potential
      if(ires.eq.1) go to 2

      if(p(1).eq.'WALLS') then
        if(n.gt.1) rmin=rlread(p(2))
        if(n.gt.2) rmax=rlread(p(3))

      elseif(p(1).eq.'SONLY') then
        bothsz=.false.
      elseif(p(1).eq.'TYPE') then
c input type of particles and lambda. There must be 2 TYPE lines in data.
        ntypes=ntypes+1
        hbs2m=rlread(p(3))
        if(hbs2m.lt.0.) then
          write (*,*)' hbs2m negative ! ',hbs2m
          stop
         endif
c effective quantumness parameter
        ehbs2m=ehbs2m+hbs2m

      elseif(p(1).eq.'GRID') then
c input radial grid
c save away grid
        do i=1,n
          grid(i)=p(i)
        enddo
        ngrid=n
c input grid parameters
        nx=intread(p(2))
c setup grid
        call setgrid(mx,nx,rv,grid(3))
        write (*,*)' grid nx=',nx,'/',mx,rv(1),rv(nx)

      elseif(p(1).eq.'SQUARER') then
c input parameters having to do with squaring
c 2. tlow is the  lowest physical temperature for which to calculate density matrix
        tlow=rlread(p(2))
        if(tlow.le.0.) then
           write (*,*) ' temperature out of range = ',tlow
           stop
        endif
        tauhigh=1./tlow
c 3.  ntemp is number of temperatures to be written out.
        ntemp=intread(p(3))
        if(ntemp.lt.1.or.ntemp.gt.mtemp)then
           write (*,*)' ntemp=',ntemp,'  out of range '
           stop
        endif
        taulow=tauhigh/2**(ntemp-1)
        write (*,*)' lowest temp =',tlow,' highest temp =',1./taulow
     +,' number of temps output=',ntemp,'/',mtemp
c 4. spatial dimensionality (1,2,3) of the interaction
        ndim=intread(p(4))
        if(ndim.le.0.or.ndim.gt.3)then
          write (*,*) 'error in ndim=',ndim
          stop
        endif
        write (*,*) ' spatial dimension=',  ndim
        if(ndim.eq.1)bothsz=.false.

c 5. norder is the highest polynomial fit to be generated (0=diagonal end-point)
        norder=intread(p(5))
        if(norder.gt.morder.or.norder.lt.0)then
          write (*,*) ' norder out of range ',norder
          stop
        endif

c 6. number of partial waves
        nl=intread(p(6))
        if(ndim.eq.1) nl=0
        if(nl.gt.ml.or.nl.lt.0) then
            write (*,*)' nl out of range =',nl
            stop
        endif
        write (*,*)' number of partial waves= ',nl,'/',ml

c 7. number of squarings from lowest temperature
c    nsquare .lt.0 is flag for primitive approximation
        nsquare=intread(p(7))
      elseif(p(1).eq.'DUMP') then
       ifdump=1
        ndump=min(n-1,5)
        if(ndump.gt.0) then
          do id=1,ndump
            idump(id)=intread(p(id+1))
          enddo
        endif

      elseif(p(1).eq.'RANK') then
c Input potential. Disregard lines between RANK and BEGIN. Subsequent data is overwritten
        call findbegn(1)
        read (1,*) (pott(i),i=1,nx)
        go to 3

      endif
      go to 1
2     write (*,*)' No potential read in'
c write zero table
      grid(2)='1 '
      write (1,*)'RANK 2 ',nx,1
      call echo(1,grid,ngrid)
      write (1,*)'LABEL 1 r'
      write (1,*)'BEGIN potential 0'
      write (1,888) (pott(i),i=1,nx)
3     continue

c check to see if all other data is present before continuing
        if(nx.eq.0)then
          write (*,*)' no grid specified  nx=',nx
          stop
        endif
        if(tlow.le.0.) then
          write (*,*)' no SQUARER line present tlow =',tlow
	  stop
        endif
	if(ntypes.eq.0) then
	  write (*,*)' ntypes=0 ehbs2m=1'
	  ehbs2m=1.d0
	elseif(ntypes.eq.1) then
	  write (*,*)' ntypes =1 DOUBLE '
	  ehbs2m=ehbs2m+ehbs2m
        endif
        if(ehbs2m.le.0)nsquare=-1
	write (*,*)' hbar^2(1/m1+1/m2)=',ehbs2m
        if(nsquare.le.0) then
            write (*,*) 'primitive approximation used ehbs2m=',ehbs2m
            norder=0
        else
            write (*,*) ' number of squarings =',nsquare
            if(nsquare.lt.ntemp) then
                write (*,*) 'nsquare.lt.ntemp'
                stop
             endif
        endif

        if(ndim.gt.1)rmin=max(rmin,0.d0)
        write (*,*)' density matrix zeroed at ',rmin,rmax
        if(rmin.ge.rv(1).or.rmax.le.rv(nx)) then
          write (*,*)' walls overlap with grid '
          stop
        endif
        write (*,*)' order of polynomial fit= ',norder,'/',morder
        if(bothsz) then
           write(*,*)' 2d polynomial fit in s and z'
        else
           write(*,*)' 1d polynomial fit in s only'
        endif
        write (*,*)' number of derivatives md=',md
        if(md.lt.1.or.md.gt.3)stop
        write (1,*)'NDERIV ',md

c****end of input section***********************************

      if(nsquare.le.0) then
c primitive approximation
          call primitive(nx,rv,pott,ntemp,taulow,ehbs2m,diag
     +    ,mx,ndim,md)

      else

c  matrix squaring
         tau=1.d0/(2.d0**nsquare*tlow)
         call initialu(nx,ml,mx2,rv,pott,tau,vint,uold,md)
         write (*,*)' Beginning matrix squaring nsquare = ',nsquare

         do isq=1,nsquare
c do a single matrix squaring
           call sqit(nx,nl,ml,mx2,ndim,tau,ehbs2m,rv,rmin,rmax,vint
     +     ,uold,unew,di,dj,md)
c double tau
           tau=tau+tau
           itemp=isq-nsquare+ntemp
           if(itemp.gt.0) then
              isd=itemp
           else
              isd=1
           endif

c compute diagonal density matrix
             call diagonal(nx,ndim,tau,ehbs2m,rv,lmax,uold,lpos,nl
     +              ,mx,diag(1,1,isd),ml,mx2,pott,md)
c write out high and low values and compute sum rule
            sa=.5d0*sang(ndim)
            do id=0,md
             sum(id)=0.d0
            enddo

            do id=1,nx
              rho=exp(-diag(id,1,isd))
              wt=sa*rv(id)**(ndim-1)

	      if(id.eq.1) then
              write (*,801) rv(id),rho,(rho*diag(id,k,isd),k=2,md)
              wt=wt*(rv(2)-rv(1))
	      elseif(id.eq.nx) then
              write (*,801) rv(id),rho-1.d0,(rho*diag(id,k,isd),k=2,md)
              wt=wt*(rv(nx)-rv(nx-1))
	      else
	      wt=(rv(id+1)-rv(id-1))*wt
	      endif
!             sum(0)=sum(0)+wt*(rho-1.d0)
              sum(0)=sum(0)+wt*(rho)
              sum(1)=sum(1)+wt*rho*pott(id)
	      if(md.ge.2)sum(2)=sum(2)+wt*rho*(diag(id,2,isd)-pott(id))
	      if(md.ge.3)sum(md)=sum(md)+wt*rho*diag(id,md,isd)
            enddo
801         format(' r ',e12.5,' rho ',3e14.7)
            write (*,802) (sum(id),id=0,md)
802         format(' cluster pf=',e12.5,' pot=',e12.5,' kes =',2e12.5)
            write (*,*) ' ke=',sum(2)/sum(0),' pe=',sum(1)/sum(0),' e='
     &  ,(sum(2)+sum(1))/sum(0)

           write (*,*) ' partial waves needed =',lmax
           if(itemp.gt.0) then 
c now do fits to polynomials
             if(norder.gt.0)call fitter(uold,diag(1,1,isd),rv,lpos
     +      ,fitdm(1,1,1,1,itemp),chis(1,0,1,itemp),norder
     +      ,ndim,ml,mx2,mx,morder2,morder,nx,rmin,rmax,tau,ehbs2m
     +      ,bothsz,md,bf(1,itemp))
!           call backflow(uold,diag(1,1,isd),rv,lpos)
c subtract potential from beta derivative
             if(md.ge.2) then
             do i=1,nx
                 diag(i,2,itemp)=diag(i,2,itemp)-pott(i)
             enddo
             endif
c readjust nl?
c           nl=min(nl,lmax+5) 
           endif
         write (*,*) '*****finished isq=',isq,' temp= ',  1./tau,'****'
	 write (*,*)' '
         enddo
      open (55,file=fname(1:ln)//'.bf',status='unknown')
         do i=1,nx
          write (55,'(8e13.5)') rv(i),(bf(i,itemp),itemp=1,ntemp)
         enddo

      endif

      
c pointers for fitting expansion
        do  k=0,norder
         if(bothsz) then
           nkt(k)=(k+3)*k/2
         else
           nkt(k)=k
         endif
        enddo
c now write out all the fitted data
      do 2000 id=1,md
       do 2000 k=0,norder
        ns=nkt(k)+1
        write (1,*)'RANK 3 ',nx,' ',ns,' ',ntemp
        grid(2)='1 '
        call echo(1,grid,ngrid)
        write (1,*)'LABEL 1 r'
        write (1,889)  3,taulow,tauhigh
889     format(' GRID ',i3,' LOG ',2e17.9)
        write (1,*)'LABEL 2 terms'
        write (1,*)'LABEL 3 time'
        write (1,*)'BEGIN density matrix table fit=',k,' derv = ',id
        do 2000 l=1,ntemp
c first diagonal
         write (1,888) (diag(i,id,l),i=1,nx)
c now off-diagonal
          do j=1,nkt(k)
          write (1,888) (fitdm(i,j,k,id,l),i=1,nx)
          enddo
888      format(5e17.9)
2000  continue
c write out rms error for the fits
      open (94,file=fname(1:ln)//'.ch',status='unknown')
      write (94,*)' RANK 4 ',nx,norder+1,md,ntemp
      write (94,*)'LABEL 1 r'
      call echo(94,grid,ngrid)
      write (94,*)'LABEL 2 order'
      write (94,*)'LABEL 3 u and du/dt'
      write (94,*)'LABEL 4 time'
      write (94,889)  4 ,taulow,tauhigh
      write (94,*)' BEGIN chi'
      write (94,888)((((chis(i,k,id,l),i=1,nx),k=0,norder),id=1,md)
     +   ,l=1,ntemp)

      if(ifdump.ne.0) then
         do id=1,md
         do i=1,nx
           write (1+id,800) rv(i),(diag(i,id,l),l=1,ntemp)
         enddo
         enddo
         do l=1,ndump
         do i=1,nx
           ind=isym(i,i)
           write (10+l,800) rv(i),(uold(idump(l),ind,id),id=1,md)
         enddo
         enddo
800      format(8e13.5)
      endif

c calculate and write out sampling tables
      call smear(ndim,ntemp,ehbs2m,taulow,diag,rv,nx,mx,md)
      close(1)
      write (*,*)' SQUARER completed successfully'
      end
