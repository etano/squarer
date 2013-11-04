       subroutine fitpnnew(v,y,t,rk,wt,nk,nf,m,mp,rad,nknots,t0,vt0,t1,
     &  vt1,coul,vmad,vol,s,rknot,delta,ndim)
! code of Natoli-Ceperley method, by Simone Chiesa, 2005
!   // ndim: spatial dimension, only 2 and 3 implemented!
!
!   // v = function in k space (INPUT)
!   // y = function in k space after break up (OUTPUT)
!   // t = coefficient of polynomials (OUTPUT)
!===== k space info =====
!   // rk = value of |k| for each shell (INPUT)
!   // wt = weight of shell (INPUT)
!   // nk = number of shell (INPUT)
!   // nf = number of shells after breakup (INPUT)
!===== splines info =====
!   // m = n of continuous derivatives [m+1 is the number of polynomial per knot] (INPUT)
!   // rad = radius of short range function (INPUT)
!   // nknots = number of knots in the interval (INPUT)
!   // t0 = constraint on the value at the origin ? (INPUT)
!   // vt0 = value at the origin (if t0=true) (INPUT)
!   // t1 = constraont on the 1st derivative at the origin ? (INPUT)
!   // vt1 = value of 1st derivative at origin (if t1=true) (INPUT)
!   // coul = want factor 1/r ? (INPUT)
!===== debug =====
!   // vmad = exact madelung constant (INPUT)

!===== INTERNAL ======
!   // delta = distance between knots
!   // s = coefficient defining polynomials
!   // mp = dimension of matrix m
!   // maxn = degree of polynomials necessary to ensure continuity of m derivatives
!   // nalphaknot = total number of polynomials

       implicit none
       integer m,nalphaknot,nknots,k,nk,i,alpha,ialpha,maxn,n,j,istart
       integer iskip,nalc,jrow,icol,nf,ld,info,lwork,ialc,mp,ndim
       double precision delta,vt0,vt1,vmad,vmad1,vmad2,chisq,vol,r,rad
       double precision s(0:mp,0:2*mp+1),ddplus(0:2*m+1)
     .                 ,ddminus(0:2*m+1)
       double precision c(0:nknots*(m+1)-1,0:nk)
       double precision v(0:nk),wt(0:nk),y(0:nk)
       double precision a(0:nknots*(m+1)-1,0:nknots*(m+1)-1)
       double precision w(0:nk),t(0:nknots*(m+1)-1),b(0:nknots*(m+1)-1)
       double precision bb(nknots*(m+1)),x(nknots*(m+1))
       double precision aa(nknots*(m+1),nknots*(m+1)),wknot(0:nknots)
       double precision rk(0:nk)
       logical coul,t0,t1
       double precision sv(nknots*(m+1)),u(nknots*(m+1),nknots*(m+1))
       double precision vt(nknots*(m+1),nknots*(m+1))
       double precision work(2000000),condnum,thrcond,wub
       parameter (thrcond=1.0d-6)  
       integer ipiv(nknots*(m+1))
       logical svd,linpack,debug
       double precision rknot(0:nknots),rpos,ww
       integer numint
       double precision pi,rcond,z(3*nknots*(m+1)-3)
       double precision dlansy,anorm
       integer z2(nknots*(m+1)-1)
      parameter  (pi=3.14159265359d0)
 
!      write (*,*) (v(k),k=0,3)
!      write (*,*) (rk(k),k=0,3)
!      write (*,*) (wt(k),k=0,3)
!      write (*,*) nk,nf,rad,nknots
!      write (99,*) t0,t1,coul,vmad,ndim
!      write (99,*) m,mp,vt0,vt1,vol
!     write (99,'(i4,3e18.7)')(i,rk(i),wt(i),v(i),i=0,99)
!     write (99,'(i4,3e18.7)')(i,rk(i),wt(i),v(i),i=nk-99,nk)


!   // LEADING DIMENSION OF aa,bb
       ld=nknots*(m+1)
       wknot(0)=0.d0
       do i=1,nknots
        wknot(i)=1.d0
       enddo

!   // DISTANCE BETWEEN KNOTS
       delta=rad/dble(nknots)

!      check if distance is too small
!
       if(rk(nk)*delta.le.20.d0*pi) then
           write(*,*)'fitpnnew: too many points',delta*rk(nk),nknots
     &,rk(nk),rad
	   stop
       end if

!   // GIVEN A DEGREE m, COMPUTE THE COEFFICIENTS DEFINING
!   // THE SPLINE FUNCTIONS
       call basis(m,mp,s)
!       write(*,*) 'COEFFICIENTS OF THE BASIS FUNCTIONS'
!       write(*,*)
!       do i=0,m
!        write(unit=*,fmt=50) (s(i,j),j=0,2*m+1)
!       enddo
!50     format(20(f14.6,1x))


!   // TOTAL NUMBER OF BASIS FUNCTIONS. 
!   // m+1 is the number of polynomial of degree 2m+1 "centered" on each knot
!   // nknots is the number of knots (from 0 to nknots-1)
       nalphaknot=nknots*(m+1)-1
!   // Degree of polynomial
       maxn=2*m+1
       write(unit=*,fmt=51) 'Using polynomials of degree',maxn,'.'
       write(unit=*,fmt=51) 
     &  '(Required for having',m,' continuous derivatives)'
       write(*,*)
 51     format(A,1x,i2,1x,A)


!   // COMPUTE THE FOURIER TRANSFORM OF SPLINES.
       do k=0,nk
        ialpha=-1
        do i=0,nknots-1
         r=i*delta
         if(ndim.eq.3) then
           call splint3D(maxn,r,delta,rk(k),vol,coul,ddplus,ddminus)
         else if(ndim.eq.2) then
           call splint2D(maxn,r,delta,rk(k),vol,coul,ddplus,ddminus)
         else
           write(6,*)'dimension not implemented in fitpnnew! ndim=',ndim
           stop
         end if
         do alpha=0,m
!      // IALPHA=(i,alpha) AND IT'S ORDERED AS (0,0) (0,1)....
!      // i INDICATES THE KNOT IN WHICH THE SPLINE IS CENTERED
!      // alpha INDICATES ONE OF THE ELEMENTARY BASIS FUNCTIONS
!      // SO ialpha=i*(m+1)+alpha (0,1,.....,nalphaknot)
          ialpha=ialpha+1
          c(ialpha,k)=0.d0
!      // COMPUTE THE F.T. of the alpha POLYNOMIAL CENTERED ON r
!      // DIFFERENT DEFINITION FOR THE KNOT r=0 (wknot=0 for r=0 otherwise =1)
          do n=0,maxn
            c(ialpha,k)=c(ialpha,k)+s(alpha,n)
     &      *(ddplus(n)+wknot(i)*ddminus(n)*(-1)**(n+alpha))
          enddo
          c(ialpha,k)=delta**alpha*c(ialpha,k)
         enddo
        enddo
       enddo
!       write(*,*) 'F.T. of all polynomials in all intervals completed'
       
!   // DEFINE THE MATRICES OF THE SYSTEM.
!   // (a(i,j) IS SYMMETRIC, I DIDN'T EXPLOIT THAT)
       do j=0,nalphaknot
        b(j)=0.d0
        do k=nf+1,nk
         b(j)=b(j)+c(j,k)*v(k)*wt(k)
        enddo
        do i=0,nalphaknot
         a(j,i)=0.d0
         do k=nf+1,nk
          a(j,i)=a(j,i)+c(j,k)*c(i,k)*wt(k)
         enddo
        enddo
       enddo

!       write(*,*) 'Matrices for the linear system generated.'

!   // SET ISTART,ISKIP,t(0) AND t(1) ACCORDING TO CONSTRAINTS
!   // ALSO CORRECT b(i) ACCORDING TO CONSTRAINS
       if (t0.and..not.t1) then
        istart=1
        iskip=0
        nalc=nalphaknot
        do j=0,nalphaknot
         b(j)=b(j)-vt0*a(j,0)
        enddo
        t(0)=vt0
        write(*,*) 'Constraint on the value active.',vt0
       elseif (t1.and..not.t0) then
        istart=0
        iskip=1
        nalc=nalphaknot
        do j=0,nalphaknot
         b(j)=b(j)-vt1*a(j,1)
        enddo
        t(1)=vt1
        write(*,*) 'Constraint on the derivative active.',vt1
       elseif (t0.and.t1) then
        istart=2
        iskip=0
        nalc=nalphaknot-1
        do j=0,nalphaknot
         b(j)=b(j)-vt0*a(j,0)
         b(j)=b(j)-vt1*a(j,1)
        enddo
        t(0)=vt0
        t(1)=vt1
        write(*,*) 'Constraints on the value and derivative active.'
     &  ,vt0,vt1
       else
        istart=0 
        iskip=-1 
        nalc=nalphaknot+1
        write(*,*) 'No constraint active.'
       endif

!   // LOAD THE MATRICES THAT ENTER THE LINEAR SYSTEM.
!   // (BETTER IF I HAD USED POINTERS)
       jrow=0
       do j=istart,nalphaknot 
        if (j.ne.iskip) then
         jrow=jrow+1
         bb(jrow)=b(j)
         icol=0
         do i=istart,nalphaknot
          if (i.ne.iskip) then
           icol=icol+1
           aa(jrow,icol)=a(j,i)
          endif
         enddo
        endif
       enddo
        write(*,*) 'Constrained matrices loaded.'

!   // SOLVE THE SYSTEM
      svd=.false.
      linpack=.true.

      if (svd) then
!   // FIRST OF ALL, DECOMPOSE aa USING SVD
       lwork=10000
       call dgesvd('A','A',nalc,nalc,aa, 
     &    ld,sv,u,ld,vt,ld,work,lwork,info)
       if (info.eq.0) then
        write(*,*) 'Singular value decomposition OK'
        write(*,*) 'Optimal  value of lwork =', int(work(1))
        write(*,*) 'Employed value of lwork =',lwork
       else
        write(*,*) 'Problems with SVD. STOP'
        stop
       endif
!   // CHECK CONDITION NUMBER
       condnum=sv(nalc)/sv(1)
       write(*,*) 'Condition number = ',condnum
       ialc=nalc
       if (condnum.lt.thrcond) then
        write(*,*) 'Condition number below threshold'
18      continue   
         ialc=ialc-1
         condnum=sv(ialc)/sv(1)
         if (condnum.gt.thrcond) goto 19
        goto 18
19      continue
        write(*,*) 'Retaining only',ialc,' sing.values out of', nalc
       endif
!   // COMPUTE THE MINIMUM LENGTH SOLUTION (NUMERICAL RECIPIES, Eq.2.9.7)
       do j=1,nalc
        x(j)=0.d0
        do i=1,ialc
         wub=0.d0
         do k=1,nalc
          wub=wub+u(k,i)*bb(k)/sv(i)
         enddo
         x(j)=x(j)+vt(i,j)*wub
        enddo
       enddo
      elseif (linpack) then
       call dpotrf('U',nalc,aa,ld,info)
!      call dpoco(aa,ld,nalc,rcond,z,info)
!      call spofa(aa,ld,nalc,info)
       if(info.ne.0) then
!          write (*,*)' inversion error dpoco info=',info
           write (*,*)' inversion error dpotrf info=',info
       endif
       anorm = dlansy('1','U',nalc,aa,ld,z)
       write(*,*) 'anorm: ',anorm
       call dpocon('U',nalc,aa,ld,anorm,rcond,z,z2,info)
       write(*,'(a26,e10.5)')'numerical inversion error ',1.d0/rcond
       if(rcond.lt.1.e-7)write(*,*)' potentially singular matrix'
       call dpotrs('U',nalc,1,aa,ld,bb,ld,info)
!      call dposl(aa,ld,nalc,bb)
!      call sposl(aa,ld,nalc,bb)
       do jrow=1,nalc
         x(jrow)=bb(jrow)
       enddo
      else
!      call dgesv(nalc,1,aa,ld,ipiv,bb,ld,info)
       do jrow=1,nalc
         x(jrow)=bb(jrow)
       enddo
      endif

!   // FILLING IN THE VALUE OF THE COEFFICIENTS
       jrow=0
       do j=istart,nalphaknot
        if (j.ne.iskip) then
         jrow=jrow+1
         t(j)=x(jrow)
        endif
       enddo
     
!   // DEBUG
       debug=.false.
       if (debug) then
        ialpha=0
        do j=0,nknots-1
         do alpha=0,m
          write(*,*) t(ialpha)
          ialpha = ialpha+1
         enddo
         write(*,*)
        enddo
       endif

!   // PLOTTING UTILITY

!   // DEFINE THE GRID
       rknot(0)=0.d0
       do i=1,nknots
        rknot(i)=rknot(i-1)+delta
       enddo

!   // COMPUTE THE VALUE OF THE FUNCTION 
!       numint=500
!       do j=1,numint-1
!        rpos=j*rad/dble(numint)
!        call computespl(rpos,ww,m,mp,maxn,nknots,s,rknot,t,delta,coul)
!        write(fmt=*,unit=55) rpos,ww*rpos,ww
!       enddo
!       write(fmt=*,unit=55) rad,0.d0

!   // DETERMINING THE FOURIER COEFFICIENT OF THE LONG RANGE PART
       do k=0,nk
        w(k)=0.d0
        do i=0,nalphaknot
         w(k)=w(k)+t(i)*c(i,k)
        enddo
        y(k)=v(k)-w(k)
       enddo
!       write(*,*) 'FOURIER TRANSFORM OF LONG RANGE DETERMINED'
 
!   // COMPUTING CHI-SQUARE
       chisq=0.d0
       do k=nf+1,nk
        chisq=chisq+y(k)**2*wt(k)
       enddo
       write(*,*) 'Chi =',sqrt(chisq)

!   // MADELUNG CONSTANT COMPUTATION 
       vmad1=t(1)
       vmad2=0.d0
       do k=0,nf
        vmad2=vmad2+wt(k)*y(k)
       enddo
       if (coul) then
        write(*,*) 'MADELUNG CONSTANT'
        write(*,*)' Computed madelung constant=',vmad1+vmad2
        write(*,*)' Gaussian charge value=',vmad
        write(*,*)' Difference ', vmad-vmad2-vmad1
        vmad=vmad-vmad2
       else
        write(*,*) 'Correction at the origin =',vmad2
       endif
       vmad=vmad1 !we need this for the non square (or cubic) case DMC
!      write (*,*) (y(k),k=1,nf)
!      write (*,*) (t(k),k=1,10)
       write (*,*) ' vmad=',vmad
       end
