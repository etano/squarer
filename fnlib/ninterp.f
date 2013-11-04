      subroutine interp(r,ix,a1,a2,a3,a4)
      include 'mach.p'
      character gtype*8
      real*8 r,a1,a2,a3,a4,dri,di,x0,rr,x,p,pm1,pm2,pp1
      integer ix,nx
      common/gcom/gtype,dri,di,x0,rr,nx
      real*8 six,si
      parameter ( six=6. , si=1./six )
c routine to produce interpolation values

      if(gtype.eq.'LINEAR')then
      x=x0+dri*r
      elseif(gtype.eq.'POWER') then
      x=x0+(dri*(r-rr))**di
      elseif(gtype.eq.'INVERSE') then
      x=x0+dri/(r-rr)
      elseif(gtype.eq.'LOG') then
      x=1+x0*log(dri*r)
      endif

c try a 3 point formula instead
      ix=x+.5
      ix=max(2,ix)
      ix=min(nx-1,ix)
      p=x-ix
      a4=0.

      if(x.lt.1) then
      a1= 1
      a2=0
      a3=0.
      elseif(x.gt.nx) then
      a1=0.
      a2=0
      a3= 1.
      else
        pm1=p-1.
        pp1=p+1.
        a1=.5*p*pm1
        a2=-pp1*pm1
        a3=.5*pp1*p
      endif
 
      return
      end
