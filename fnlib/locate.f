      subroutine locate(xx,n,x,j)
c
c NR binary search routine
      integer j,n
      real*8 x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
 10   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
            jl=jm
         else
            ju=jm
         endif
      goto 10
      endif
      j=jl
      return
      end
