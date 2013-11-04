      subroutine breakup(p,n) ! does breakup of long range potential
      implicit none


      integer mnkv,mnsh,mdim,mnshex
      parameter (mnkv=100000,mnsh=1000,mdim=3,mnshex=8600)

      integer ndim,nlamb,l0,intread,n,nshlls,kmult(0:mnsh),nvects 
     &,nshex,k,l
      character p(n)*(*)
      real*8 ell(0:3),tpiell(3),cutr,eps,dstlayer,rslayer,pi,vol,rs
     &,rlread,cutk,rknorm(0:mnshex),rkcomp(mdim,mnkv),sang,sangle
     &,vkbare(0:mnshex),xpon,epsil(0:mnshex),gm,vt0,vt1,vmad
     &,vlr(0:mnshex),vzero,wtk(0:mnshex),cutk_small
      include "cbreak.cm"
      logical coul,t0,t1

      cutr=rlread(p(2))
      cutk=rlread(p(3))
      eps=rlread(p(4))
      dstlayer=0.d0

!unpack some parameters
      if(p(1).eq.'SCCLOPT') then
       rs=rlread(p(5))
       ndim=3  ! three dimensions only
       l0=5

      elseif(p(1).eq.'COULOPT') then
        ndim=intread(p(5))
        if(ndim.lt.2.or.ndim.gt.3)stop
         if(n.lt.5) then
          write (*,*)' too few parameters for COULOPT',n
          stop
         endif
         dstlayer=0.d0
         if(n.gt.5) dstlayer=rlread(p(6))
         if(n.gt.6) rslayer=rlread(p(7))
         if(rslayer.gt.0.d0) dstlayer=dstlayer/rslayer
c default with 7 parameters is cubic box and ell=2*cutr. Otherwise read in ell.
         l0=7
      endif
 
!setup box
      cutk_small=3*cutk
      pi=3.1415 92653 58979d0
      vol=1.d0
      ell(0)=2.d0*cutr
      do l=1,ndim
        if(l+l0.le.n) then
         ell(l)=rlread(p(l+l0))
        else
         ell(l)=ell(l-1)
        endif
        vol=vol*ell(l)
        if(1.9999*cutr.gt.ell(l))stop 'ipot: 1.9999*cutr.gt.ell(l)'
        tpiell(l)=2.d0*pi/ell(l)
        cutk_small= max(cutk_small,tpiell(l)*5.)
       enddo

!generate set of k-vectors
      call shells(ndim,tpiell,cutk_small,nshlls,rkcomp,rknorm(1),kmult
     +     ,nvects,mnkv,mnsh,mdim)
      nlamb=0 !find number of k-vectors less than cutk
      do k=1,nshlls
       if(rknorm(k).lt.cutk) nlamb=k
      enddo
      write (*,*) ' cutk =',cutk,' nlamb ',nlamb
! set up flags for fitpnnew
      nknots=20  ! number of spacings in spline table
      nshex=50.*pi*nknots/(cutr*rknorm(1)) ! set limit for nshex to avoid instability in fitpnnew
      if(nshex.gt.mnshex) then
            write (*,*)'increase mnshex in breakup to be >',nshex 
     &     ,' current value=',mnshex
            stop
      endif
!     nshex=min(mnshex,nshex)
      call fillk(rknorm,wtk,kmult,nshlls,nshex,ndim,vol) ! go to large k on a grid
      sangle=sang(ndim)

!generate k-space potential
      if(p(1).eq.'COULOPT') then
        vkbare(0)=0.d0 ! uniform background
        do k=1,nshex
           vkbare(k)=(eps*sangle)/(vol*rknorm(k)**(ndim-1))
!    &     * exp(-dstlayer*rknorm(k) )  ! screening
        enddo 
c Compute madelung constant.
      xpon=1.d0
      call mdlng(ndim,xpon,ell(1),eps,vmad)
      write (*,*)' vmad ',vmad
      vmad=0.d0
      if(ndim.eq.2) write (*,*) 'ideal square madelung constant',
     &-3.90026492d0*eps/ell(1)
      if(ndim.eq.3) write (*,*) 'ideal cube madelung constant',
     & -2.837297479d0*eps/ell(1)

      elseif(p(1).eq.'SCCLOPT') then
           do k=0,nshex! change units to rs
            rknorm(k) = rknorm(k)/rs
           enddo
         call epsk1(rs,epsil(1),rknorm(1),nshex,gm) ! evaluate epsilon
          do k=0,nshex
            rknorm(k) = rknorm(k)*rs !change back
           enddo
        do k=1,nshex
          vkbare(k)=(eps*sangle)/(vol*rknorm(k)**(ndim-1)*epsil(k))
        enddo 
         write (*,*)'gamma0=',gm
c  vkbare(0)=[(k_f/k_tf)^2/rs-gamma0]/k_f^2 is the k=0 limit of the potential
         vkbare(0)=((3*pi**2/16)**(2./3.)/rs-gm)/(1.919158293/rs)**2
         vkbare(0)=(eps*sangle)*vkbare(0)
         write (*,*) 'limiting value of the potential for small k:'
     &                ,vkbare(0)
         vkbare(0)=vkbare(0)/vol
c give the interaction energy with images:
        vmad=vzero(rs)

      endif
 

      m=2
      maxn=2*m+1
      if(dstlayer.gt.0.d0) then!for offset no cusp condition
          coul=.false.
          t0=.false.
          t1=.true.
          vt1=0.d0
      else
          coul=.true.
          t0=.true.   ! fixes const term for r=0 with value vt0
          vt0=eps
          t1=.false.
          vt1=0.d0
      endif

      write (*,*) 'small k=',cutk_small,' large k=',rknorm(nshex)
     & ,' knots=',nknots,'nshex=',nshex,' ndim=',ndim,'nlamb=',nlamb
     &,' m=',m,'maxm=',maxm,'cutr=',cutr
      write (*,*)  't0=',t0,vt0
      write (*,*)  't1=',t1,vt1
      write (*,*) 'coul=',coul
      write (*,*) ' vol=',vol

      call fitpnnew(vkbare,vlr,b,rknorm,wtk,nshex,nlamb,m,maxm
     &   ,cutr,nknots,t0,vt0,t1,vt1,coul,vmad,vol,sm,rknot,delta,ndim)

      write (3,301).5d0*vmad ! vimage=1/2 limit (r=0) [V_image(r) -V_bare(r)]
301   format('VIMAGE ',e13.5)

! write out k-space potential. Unit 4 should have been opened already.
      write (4,*)'RANK 1 ',nlamb+1
      write (4,*)'BEGIN k-space potential'
      write (4,*) (vlr(k),k=0,nlamb)
      close(4)
      do k=0,nshlls
       write (8,*) rknorm(k),vlr(k)*vol
      enddo
      write (*,*)' vol =',vol,'nshlls=',nshlls
      close (8)
      end
