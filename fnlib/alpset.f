        subroutine alpset(a,rc,itol,alpha)
        include 'mach.p'
        real*8 a,rc,alp,dalp,pi,spi,al2,gr,gr0,phia,phirc
        real*8 dd,t1,t2,al2a,aa,alpha
        integer isign,isign0,icon,itol

        pi = 4.*atan(1.)
        spi = 2./sqrt(pi)
        icon = 1
        isign = 1
        isign0 = isign
        t1 = 10.**(-itol)
        t2 = t1*0.9
        alp = 1.
        dalp = 0.5
c       print *,t1
1       if (isign.ne.isign0) then
         icon = icon + 1
         dalp = isign*abs(dalp)/icon
        endif
        alp = alp + dalp
c       print *,alp
        al2a = (alp*a)**2
        call gammi(gr,.5,al2a,gr0)
        phia = gr/a/gr0
        al2 = (alp*rc)**2
        call gammi(gr,.5,al2,gr0)
        phirc = gr/rc/gr0
c       fa = - phia/a - alp*exp(-al2a)/a*spi
c       frc = - phirc/rc - alp*exp(-al2)/rc*spi
c       print *,phia/a,fa,phirc/rc,frc
c       if (abs(fa).le.1.e-20) then
        if (abs(phia).le.1.e-20) then
         alp = alp*0.1
         goto 1
        end if
c       dd = frc/fa
        dd = phirc/phia
        if (dd.lt.1.e-50) then
c       alp is two large to converge
         alp = alp/5.
         go to 1
        end if
        isign0 = isign
        if (dd.gt.t1) then
         isign = 1
         goto 1
        else if (dd.le.t2) then
         isign = -1
         go to 1
        else
         alpha = alp
         write (*,*) 'subroutine alpset'
         write (*,*) 'alpha = ',alp
c        write (*,*) 'frc/fa =',dd
         write (*,*) 'phirc/phia = ',phirc/phia
        endif
        aa = 2.*log(1.e-3*rc/a)/(a*a-rc*rc)
c       print *,aa
        return
        end
