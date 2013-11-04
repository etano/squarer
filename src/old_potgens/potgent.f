
      program potgen
      implicit real*8 (a-h,o-z)
c************************************************************************
c	May 12, 1997
c	Include Ko(r) potential -- relevant for vortices
c************************************************************************
      parameter (mn=10,mx=1000)
      dimension pott(mx),rv(mx),diff(mx)
      character fname*14,p(mn)*28,eunit*3,lunit*3
     +,pot(mn)*28,grid(mn)*28

      write (*,*)' input prefix of files '
      read (*,99) fname
99    format(a14)
      ln=index(fname,' ')-1
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
        np=n-1

      endif
      go to 1

2     continue
c check if potential and grid have been specified.
      if(ifpair.ne.1)stop
      if(nx.le.0)stop
c setup potential
      do i=1,nx
      call ipot(pot,np,pott(i),rv(i),eunit,lunit)
      enddo

c now check grid
      npts=10
      do i=1,nx-1
       diff(i)=0.
       do k=1,npts
         x=rv(i)+(rv(i+1)-rv(i))*k/(npts+1.d0)
         call ipot(pot,np,potx,x,eunit,lunit)
         call interp(x,ix,a1,a2,a3,a4)
         dp=a1*pott(ix-1)+a2*pott(ix)+a3*pott(ix+1)+a4*pott(ix+2)-potx
         diff(i)=max(diff(i),abs(dp))
       enddo
         diff(nx)=diff(nx-1)
      enddo

c zero potential at cutr
      cutr=rlread(pot(2))
       call ipot(pot,np,pot0,cutr,eunit,lunit)
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
969   format(5e17.9)
      close(3)

      open (9,file=fname(1:ln)//'.dg',status='unknown')
c     write (9,*)'RANK 1 ',nx
c     call echo(9,grid,ng)
c     write (9,*)'LABEL 1 r'
c     write (9,*)'BEGIN grid error'
c     write (9,969) (diff(i),i=1,nx)
      do i=1,nx
      write (9,969) rv(i), diff(i),pott(i)
      enddo
      close(9)
      end

      subroutine ipot(p,n,pott,rv,eunit,lunit)


c interprets potential parameter and sets up table.
      character p(n)*(*),eunit*(*),lunit*(*),cite*80
      real*8 pot(3)
      save icall,cutr
      data icall/0/
      icall=icall+1
      if(icall.eq.1)ntail=0

      if(p(1).eq.'BESSK0') then
         pott=bessk0(rv)
      elseif (p(1).eq.'GAUSS') then
	 do i = 2, n
	    pot(i-1) = rlread(p(i))
	 enddo
	 np = n-1
	 pott = gauss(rv,pot,np)
      else
         write (*,*)' this potential not defined ',p(1),n
         stop
      endif
      end

      FUNCTION bessi0(x)
      DOUBLE PRECISION  bessi0,x
      DOUBLE PRECISION  ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     *1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     *0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,
     *-0.1647633d-1,0.392377d-2/
      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        bessi0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *(q7+y*(q8+y*q9))))))))
      endif
      return
      END
c**********************************************************************
      FUNCTION bessk0(x)
      DOUBLE PRECISION bessk0,x
CU    USES bessi0
      DOUBLE PRECISION bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0d0) then
        y=x*x/4.d0
        bessk0=(-log(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
      else
        y=(2.d0/x)
        bessk0=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END

	FUNCTION gauss(x,pot,np)
	implicit real*8(a-h,o-z)
	integer np
	double precision gauss, pi,pot(np)
	sigma = pot(2)
	v0 = pot(3)
	s = 2.d0*sigma*sigma
	gauss = v0*dexp(-x*x/s)
	end
