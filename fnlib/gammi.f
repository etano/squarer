      subroutine gammi(g,a,x,gm)
      include 'mach.p'
c incomplete gamma function=g(a,x)  gm=g(a,0)
      real*8 g,a,x,gm,xcut,y,xa,gami,fn,term,gamma
      integer ns,nf,n,m
      data xcut,ns,nf/1.2,15,20/
      save xcut,ns,nf
      y=abs(x)
      xa=y**a
      gm=gamma(a)
      if(y.gt.xcut) go to 2
      gami=gm-xa/a
      fn=1.
      term=-xa
      do 3 n=1,ns
      term=-term*y/n
3     gami=gami+term/(a+n)
       g=gami
      return
2     term=2*y
      m=nf
      do 4 n=1,nf
      term=y+(m-a)*term/(term+m)
4     m=m-1
      gami=exp(-y)*xa/term
       g=gami
      return
      end
