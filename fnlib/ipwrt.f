      subroutine ipwrt(qid,ndim,nparts,nslices,el2,ntf,ncf,hbf,tau,mcp)
      include 'mach.p'

      integer ndim,nparts,nslices,ntf,ncf,mcp,ln,nwp,l,i
      real*8 el2,hbf,tau

      dimension el2(ndim),ncf(ntf),hbf(ntf)
      character qid*(*),filej*32,qtime*26
c initializes file 90 by opening it and writing a header consisting of
c  qid=run id, (type character)
c  ndim= spatial dimensionality (probably 1,2 or3)
c  nparts= # of particles
c  nslices = number of time slices
c  el2= 1/2 box size in each direction
c  ntf=number of different types of particles
c  ncf(i)=#number of particles of this type
c  hbf=hbar**2/2*mass for each  type of particle
c tau=time slice
       ln=index(qid,' ')-1
       filej=qid(1:ln)//'.pc'
       open(unit=90,file=filej,form='formatted')
!      open(unit=90,file=filej,form='unformatted')
       call timedate(qtime)
c  nwp will be number of packed words/record written in wpack
      nwp=1+(ndim*nparts*nslices-1)/mcp
!     write (90) qid,qtime,ndim,nparts,nslices,(el2(l),l=1,ndim)
!    +,ntf,(ncf(i),i=1,ntf),(hbf(i),i=1,ntf),tau,nwp
      write (90,'(i1,i6,3e15.6)') ndim,nparts,(el2(l),l=1,ndim)
      end
