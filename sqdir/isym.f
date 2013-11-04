      function isym(i,j)
      implicit none
c routine to find index in packed table
      integer ipoint,isym,i,j,i1,i2
      common/ipp/ipoint(1000)
      i1=amin0(i,j)
      i2=i+j-i1
      isym=i1+ipoint(i2)
      return
      end
