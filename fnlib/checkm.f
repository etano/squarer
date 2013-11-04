      subroutine checkm(n2,a,b,it,ich)
      include 'mach.p'
      integer n2,it,ich,i,j
      real*8 a,b,bm
      dimension a(n2),b(n2)
      bm=0.d0
      do 1 i=1,n2
1     bm=max(bm,abs(b(i)-a(i)))
      if(bm.gt.1.e-3) then
        write (6,*) ' checkm ',bm,n2
        write (6,*) ' it ',it,' ich ',ich
        write (6,*) ' new matrix '
        write (6,4) (a(j),j=1,n2)
        write (6,*) ' old matrix '
        write (6,4) (b(j),j=1,n2)
4       format(5e12.5)
        write (6,*) 'new-old'
        write (6,4) (a(j)-b(j),j=1,n2)
        stop
      endif
      end

      subroutine checkm4(n2,a,b,it,ich)
      include 'mach.p'
      integer n2,it,ich,i,j
      real*4 a,b,bm
      dimension a(n2),b(n2)
      bm=0.d0
      do 1 i=1,n2
1     bm=max(bm,abs(b(i)-a(i)))
      if(bm.gt.1.e-3) then
        write (6,*) ' checkm ',bm,n2
        write (6,*) ' it ',it,' ich ',ich
        write (6,*) ' new matrix '
        write (6,4) (a(j),j=1,n2)
        write (6,*) ' old matrix '
        write (6,4) (b(j),j=1,n2)
4       format(5e12.5)
        write (6,*) 'new-old'
        write (6,4) (a(j)-b(j),j=1,n2)
        stop
      endif
      end
