      subroutine dtabl1(m,i,p,t,f,fp,n,dri)
      include 'mach.p'
c finds f and dlogf/dr from the table in t
c
      integer  m,i,n,k
      real*8 p,t, f,fp,dri
      dimension i(m),p(m),t(n,4),f(m),fp(m)
      do 1 k=1,m
      f(k)=t(i(k),1)+p(k)*(t(i(k),2)+p(k)*(t(i(k),3)+p(k)*t(i(k),4)))
      fp(k)=dri*(t(i(k),2)+p(k)*(2.*t(i(k),3)+3.*t(i(k),4)*p(k)))/f(k)
1     continue
      return
      end
