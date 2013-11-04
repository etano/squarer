      program potgen_lr 
!version of potgen only for long range potentials 5.23.07
      implicit none
      integer mn,mx
      parameter (mn=20,mx=2100)
      real*8 pott(mx),rv(mx),diff(mx),cutr,pot0,a1,a2,a3,a4,x,dp
     &,rlread,potx
      character fname*28,p(mn)*28,eunit*3,lunit*3
     +,pot(mn)*28,grid(mn)*28
       integer ntypes,nx,ifpair,j,ipickoff,n,ng,ix,i,k,npts
     &,intread,ln,n1

      character *100 BUFFER
      call getarg(1,BUFFER)
      read(BUFFER,*) fname
      write (*,*)' input prefix of files '
c      read (*,99) fname
99    format(a20)
      ln=index(fname,' ')-1
      write (*,*) fname
      open(1,file=fname(1:ln)//'.in')
c output file
      open (3,file=fname(1:ln)//'.dm',status='unknown')

      ntypes=0
      nx=0
      ifpair=0
c loop thru reading the input
1     continue
      j=ipickoff(1,p,n,mn)
      if(j.eq.1) go to 2
c repeat input on output
      call echo(3,p,n)

c input energy units and length units
      if(p(1).eq.'UNITS') then
      eunit=p(2)
      lunit=p(3)
 
      elseif(p(1).eq.'GRID')then
c input grid parameters
        do i=1,n
        grid(i)=p(i)
        enddo
        ng=n
        nx=intread(p(2))
        call setgrid(mx,nx,rv,grid(3))

      elseif(p(1).eq.'POT') then
c  input type of potential
        ifpair=ifpair+1
c store away the potential specification
        do i=2,n
         pot(i-1)=p(i)
        enddo
        n1=n-1
        open (4,file=fname(1:ln)//'.vk',status='unknown')
        open (8,file=fname(1:ln)//'.yk',status='unknown')

      endif
      go to 1

2     continue
c check if potential and grid have been specified.
      if(ifpair.ne.1)stop
      if(nx.le.0)stop

c do breakup and put results in common cbreak.cm
      call breakup(pot,n1)

      do i=1,nx ! evaluate potential
          call getpot(pott(i),rv(i))
      enddo
!     write (*,*) 'exit getpot'

c now check grid
      npts=10
      do i=1,nx-1
       diff(i)=0.d0
       do k=1,npts
         x=rv(i)+(rv(i+1)-rv(i))*k/(npts+1.d0)
         call getpot(potx,x)   
         call interp(x,ix,a1,a2,a3,a4)
         dp=a1*pott(ix-1)+a2*pott(ix)+a3*pott(ix+1)+a4*pott(ix+2)-potx
         diff(i)=max(diff(i),abs(dp))
       enddo
         diff(nx)=diff(nx-1)
      enddo

c zero potential at cutr
      cutr=rlread(pot(2))
      call getpot(pot0,cutr)
      write (3,33)pot0
33    format(' POTTAIL ',e14.5,i3,5(f7.2,e13.5))
       do i=1,nx
       if(rv(i).le.cutr) then
         pott(i)=pott(i)-pot0
       else
         pott(i)=0
       endif
      enddo
 
      grid(2)='1 '
      write (3,*)'RANK 2 ',nx,1
      call echo(3,grid,ng)
      write (3,*)'LABEL 1 r'
      write (3,*)'BEGIN potential 0'
      write (3,969) (pott(i),i=1,nx)
969   format(5e18.10)
      close(3)

      open (9,file=fname(1:ln)//'.dg',status='unknown')
      do i=1,nx
      write (9,969) rv(i),pott(i),diff(i)
      enddo
      close(9)
      end

      subroutine getpot(pott,rv) ! now evaluate potential
      implicit none
      real*8 rv,pott,dpot,ddpot,d3pot
      include "cbreak.cm"

!     write (*,*) 'rv',rv,m,maxm,maxn,nknots,delta
      call computespl2(rv,pott,dpot,ddpot,d3pot,m
     .           ,maxm,maxn,nknots,sm,rknot,b,delta)
      pott=pott/rv

      end
