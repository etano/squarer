      subroutine ggrid(r,style,n,r0,r1)
      include 'mach.p'
c generates a mesh of "style"; n=number of points
c r0 is first point r1 last point

      real*8 r,r0,r1,dr,rr
      integer n,i

      dimension r(n)
      character style*(*)

      if(style.eq.'LINEAR') then
        dr=(r1-r0)/(n-1)
        do 10 i=1,n
10      r(i)=r0+(i-1)*dr
      elseif(style.eq.'LOG') then
        dr=(r1/r0)**(1./(n-1.))
        rr=r0
        do 20 i=1,n
        r(i)=rr
20      rr=rr*dr
      elseif(style.eq.'POLY') then

      else
c currently undefined grid style
        write(*,*) 'currently undefined grid style in ggrid: ',style
        stop
      endif
      return
      end
