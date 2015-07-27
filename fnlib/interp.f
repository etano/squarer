      subroutine interp(r,ix,a1,a2,a3,a4)
      include 'mach.p'
      character gtype*8
      real*8 r,a1,a2,a3,a4,dri,di,x0,rr,x,p,pm1,pm2,pp1
      integer ix,nx,i
      common/gcom/gtype,dri,di,x0,rr,nx,gdri,gx0,gnr !this common communicates with interp
      real*8 six,si,sumk
      parameter ( six=6.d0 , si=1.d0/six )
c routine to produce interpolation values

      if(gtype.eq.'LINEAR')then
      x=x0+dri*r
      elseif(gtype.eq.'POWER') then
      x=x0+(dri*(r-rr))**di
      elseif(gtype.eq.'INVERSE') then
      x=x0+dri/(r-rr)
      elseif(gtype.eq.'LOG') then
      x=1+x0*log(dri*r)
      elseif(gtype.eq.'LOGLIN') then
        if (r.lt.gnr) then
          x=1+gx0*log(gdri*r)
        else
          x=x0+dri*r
        endif
      endif

      ix=x
      ix=max(2,ix)
      ix=min(nx-2,ix)
      p=x-ix

      if(x.le.nx.and.x.ge.1)then
        pm1=p-1.
        pm2=p-2.
        pp1=p+1.
        a1=-si*p*pm1*pm2
        a2=.5*pm1*pp1*pm2
        a3=-.5*p*pp1*pm2
        a4=si*p*pp1*pm1
      else
        a1=(2-p)/3.
        a2=0.
        a3=0.
        a4=(p+1)/3.
      endif
 
      return
      end
