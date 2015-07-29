      subroutine setgrid(mx,nxx,r,p)
c returns a grid in r and puts constants in common
      include 'mach.p'
      real*8 r,dri,di,x0,rr,r1,rn,dr,d,gg,f1,fn,rlread,y,gdr
      integer mx,nxx,nx,i
      dimension r(mx)
      character gtype*8, p(5)*(*)
      common/gcom/gtype,dri,di,x0,rr,nx,gdri,gx0,gnr !this common communicates with interp
      gtype=p(1)
      nx=nxx
      if(nx.le.0.or.nx.gt.mx)stop
      r1=rlread(p(2))
      rn=rlread(p(3))
      if(nx.eq.1.and.r1.eq.rn) then
        r(1)=r1
        return
      endif
      if(r1.ge.rn)stop
      write (*,*) r1,rn

      if(gtype.eq.'LINEAR')then
        dri=(nx-1.)/(rn-r1)
        x0=1-r1*dri
        dr=1./dri
        do i=1,nx
        r(i)=r1+dr*(i-1)
        enddo
       write (*,*) gtype,dri,di,x0,rr,nx
 
      elseif(gtype.eq.'LOG') then
        if(rn*r1.le.0.)stop
        x0=(nx-1)/log(rn/r1)
        dri=1./r1
        dr=(rn/r1)**(1./(nx-1.))
        do i=1,nx
        r(i)=r1*dr**(i-1)
        enddo
      elseif(gtype.eq.'INVERSE') then
       y=rlread(p(4))
       if(y.eq.1.or.y.le.0)stop
       x0=(nx-y)/(1-y)
       dri=-(rn-r1)*(nx-x0)*(1-x0)/(nx-1.)
       if(dri.ge.0.)stop
       rr=r1-dri/(1-x0)
       do i=1,nx
       r(i)=rr+dri/(i-x0)
       enddo
 
      elseif(gtype.eq.'POWER') then
        x0=rlread(p(4))
        d=rlread(p(5))
        if(d.eq.0.or.x0.ge.1.)stop
        gg=(r1/rn)**(1./d)
        x0=(nx*gg-1)/(gg-1)
        f1=(1-x0)**d
        fn=(nx-x0)**d
        di=1./d
        dri=(fn-f1)/(rn-r1)
        rr=(r1*fn-rn*f1)/(fn-f1)
        dr=1./dri
        do i=1,nx
        r(i)=rr+dr*(i-x0)**d
        enddo

      elseif(gtype.eq.'LOGLIN') then
        gnr=rlread(p(4))
        gx0=((nx/2.)-1.)/log(gnr/r1)
        gdri=1./r1
        gdr=(gnr/r1)**(1./((nx/2.)-1.))
        dri=(((nx/2.)+1.)-1.)/(rn-gnr)
        x0=((nx/2.))-gnr*dri
        dr=1./dri
        do i=1,nx
          if(i>nx/2.) then
            r(i)=gnr+dr*(i-(nx/2.))
          else
            r(i)=r1*gdr**(i-1)
          endif
        enddo

      else
        write (*,*)' grid type undefined '
        write (*,*) gtype
        stop
      endif
      write (*,*) nxx,r1,rn

      return
      end
