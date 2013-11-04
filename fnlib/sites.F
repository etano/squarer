      subroutine sites(ir,x,natoms,ndim,nxtalp,vol,ncell,ell,rnn,mdim)
      implicit none
#define SIMPLE_SPRNG    
#include "sprng_f.h"

      integer ir,natoms,ndim,nxtalp,ncell(3),ii,l,nxtal,nvac
      integer npuc,nc,nvact,is,k,it,npts,icount(3),i,j,mdim
      real*8 vol,rnn,small,ro,d,factor,ds,qq,alpha
     & ,xcell,a,ellsq,rn2,rsq,xt,b(3)
     & ,x(mdim,natoms),ell(mdim),q(3,16)
      data small/1.0d-5/
      save small
c******************************************************
c  sites computes natoms crystal sites and puts them in x
c  ndim is the spatial dimensionality
c  nxtal is the crystal type see 1,2,3,4 below
c      if zero will chose crystal type to minimize number of vacancies
c  ro is the natoms per unit volume-used to compute ell
c  ncell(3) are the number of unit cells in each direction
c    if product(ncell)*npuc.lt.natoms ncell is increased so that
c    there at least as many lattice sites as particles and for hcp
c    the box is roughly cubic
c ell(3) --computed--is the size of the simulation box =ncell*cell size
c*********************************************************************
      ro=natoms/vol
      if(ir.ne.0)write (6,15) natoms,ndim,ro
15    format(' computing',i4,' lattice sites  dimensionality',i2
     +,' number density ',e12.5)
      do  ii=1,8
      do  l=1,3
         q(l,ii)=0.0d0
      enddo
      enddo
      do  l=1,ndim
         b(l)=1.0d0
      enddo
      nxtal=iabs(nxtalp)
      if(nxtal.eq.0) then
c determine lattice type by minimizing 2**l*ncell**ndim-natoms
      nvac=natoms
      do  l=1,ndim
      npuc=2**(l-1)
      nc=(float(natoms)/float(npuc))**(1./ndim)+1.0-small
      nvact=npuc*nc**ndim-natoms
      if(nvact.lt.nvac) then 
         nvac=nvact
         nxtal=npuc
      endif
      enddo
      if(ir.ne.0)write (6,461)
461   format('  crystal type chosen by default to minimize vacancies ')
      endif
      npuc=nxtal

       if(nxtal.eq.1) then
      if(ir.ne.0)write (6,*)'  simple cubic lattice'

      elseif(nxtal.eq.2) then
         if(ir.ne.0)write (6,*)'  body centered cubic lattice'
         d=.5d0
         do  l=1,ndim
           q(l,2)=d
         enddo

      elseif(nxtal.eq.3.or.nxtal.eq.9) then
      if(ir.ne.0)write (6,*)'   hexagonal close packed lattice'
      npuc=4
      if(ndim.eq.2) npuc=2
      if(ndim.eq.2) b(1)=3.d0**(-0.25d0)
      if(ndim.eq.3) b(1)=1.d0/sqrt(2.d0)
      b(2)=b(1)

!     b(2)=b(1)*dsqrt(3.d0) ! bravais vector along x (traditional)
      b(1)=b(2)*dsqrt(3.d0) !  bravais vector along y

      b(3)=b(2)*sqrt(2.d0)/1.5d0
      q(1,2)=0.5d0
      q(2,2)=0.5d0
      q(1,3)=0.5d0
      q(2,3)=5.d0/6.d0
      q(3,3)=0.5d0
      q(1,4)=1.0d0
      q(2,4)=1.d0/3.d0
      q(3,4)=0.5d0

      if(nxtal.eq.3) go to 10
      if(ir.ne.0)write(6,*)'mhcp lattice'
      npuc=8
c scale by the c/a ratio
      factor=(x(2,1)*sqrt(.375))**(-1./3.)
      b(1)=b(1)*factor
      b(2)=b(2)*factor
      b(3)=b(3)/factor**2
      ds=x(1,1)*ro**(1./3.)/(4.*b(3))
      do  is=1,4
      do  l=1,3
      qq=q(l,is)
      if(l.le.2) then
        d=0.
      else
        d=ds
      endif
       q(l,is)=qq-d
       q(l,is+4)=qq+d
      enddo
      enddo
c b is the size of the unit cell--volume of unit cell is one
c q(npuc,ndim)*b(ndim) are the vector dispacements of sites withincell

      elseif(nxtal.ge.4.and.nxtal.le.7) then
      if(ir.ne.0)write (6,*)'  face-centered cubic lattice'
      do  ii=1,3
      do  l=1,3
      if(ii.ne.l) q(l,ii+1)=.5
      enddo
      enddo
      if(nxtal.eq.4) go to 10
      npuc=8
      if(nxtal.eq.6) go to 810
      if(nxtal.eq.5) go to 801

c now add alpha nitrogen displacements
      if(ir.ne.0)write (6,*)' alpha nitrogen fcc lattice'
      ds=x(1,1)*ro**(1./3.)/sqrt(48.0)
      do  l=1,3
      do  is=1,4
        d=-ds
        if(q(l,is).gt..2) d=ds
        qq=q(l,is)
        q(l,is)=qq+d
        q(l,is+4)=qq-d
      enddo
      enddo
      go to 10

801   if(ir.ne.0) write(6,*)' diamond lattice'
      ds=.125
      do  l=1,3
      do  is=1,4
      qq=q(l,is)
      q(l,is)=qq-ds
      q(l,is+4)=qq+ds
      enddo
      enddo
      go to 10

810   if(ir.ne.0) write (6,*)' diamond bond lattice'
      npuc=16
      k=4
      ds=.25
      do  is=1,4
      do  it=1,3
         k=k+1
         do  l=1,3
         q(l,k)=q(l,is)
         if(l.ne.it)q(l,k)=q(l,k)+ds
         enddo
       enddo
       enddo


      elseif(nxtal.eq.8) then
      if(ir.ne.0)write(6,*)'simple hexagonal lattice'
      npuc=  2
      alpha=1./(x(2,1)*sqrt(3.))**(1./3.)
      b(1)=alpha
      b(2)=sqrt(3.)*alpha
      b(3)=x(2,1)*alpha
      do  l=1,2
        q(l,2)=.5
      enddo

      endif
10    continue
      if(ir.ne.0)write (6,505) ((q(l,is),l=1,ndim),is=1,npuc)
505   format(' q displs ',3f10.5)
      npts=1
      do  l=1,ndim
         npts=npts*ncell(l)
      enddo
      if(npts*npuc.ge.natoms) go to 30
c recalculate ncell since there are too few lattice points
      xcell=(float(natoms)/float(npuc))**(1.d0/float(ndim))
      npts=1
      do  l=1,ndim
        ncell(l)=(xcell/b(l))+1.0-small
        npts=npts*ncell(l)
      enddo
30    a=(vol/npts)**(1./ndim)
      do  l=1,ndim
        icount(l)=0
        ell(l)=ncell(l)*a*b(l)
      enddo
      nvac=npts*npuc-natoms
      if(ir.ne.0)write (6,255) npuc,nvac,(ell(l),l=1,ndim)
255   format(' npuc ',i2,' nvacancies ',i5,' box size',3e12.5)
      if(ir.ne.0)write (6,256) (ncell(l),l=1,ndim)
256   format(' number of cells in each direction',3i5)
c put in nvac vacancies by flags in x array
      do  i=1,natoms
         x(1,i)=0.d0
      enddo
      if (nvac.eq.0) go to 9
c compute index  for skipping
      do  i=1,nvac
        j=(natoms-1)*sprng()+2.
c note that the first site will always be filled
        x(1,j)=x(1,j)+1.
      enddo
9     continue
      i=1
c  the particles are confined abs(x(l)).le.el2(l)
      ellsq=ell(1)**2
      rnn=ellsq
      rn2=ellsq
c   loop over the different points inthe unit cell
      do  is=1,npuc
c loop over all the unit cells
      do  j=1,npts
c skip over lattice site if vacant
      if(x(1,i).lt.small) go to 8
      x(1,i)=x(1,i)-1.d0
      if(ir.ne.0)write (*,*)' vacancy ',is,' cell ',(icount(l),l=1,ndim)
      go to 300
8     rsq=0.0d0
      do  l=1,ndim
        xt=a*b(l)*(icount(l)+q(l,is))
        if(xt.ge.0.5*ell(l)) xt=xt-ell(l)
        if(xt.lt.-.5*ell(l)) xt=xt+ell(l)
        x(l,i)=xt
        rsq=rsq+(xt-x(l,1))**2
      enddo
      i=i+1
c find the nearest and next nearest distance from the origin
c exclude the origin
      if(rsq.lt.small*rnn) go to 300
c is this point greater than previously found 2nd minima
      if(rsq-rn2.gt.-small*rnn) go to 300
      if(abs(rsq-rnn).lt.rnn*small) go to 300
c rnn and rsq must be first and second smallest distances.
      rn2=rnn
      if(rsq.gt.rn2) rn2=rsq
      if(rsq.lt.rnn) rnn=rsq
300    continue
      do  l=1,ndim
        icount(l)=icount(l)+1
        if(icount(l).lt.ncell(l)) go to 450
        icount(l)=0
      enddo
450   continue
      enddo
      enddo
      rnn=sqrt(rnn)
      rn2=sqrt(rn2)
      if(ir.ne.0)then
       write (6,100) rnn,rn2
100   format(' nearest and next nearest neighbor distance  ',2e12.5)
c      do 110 i=1,natoms
c110    write (6,111) i,(x(l,i),l=1,ndim)
c111    format(' lattice sites ',i5,3e12.5)
       endif
      end
