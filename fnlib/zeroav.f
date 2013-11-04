 
      subroutine zeroav
      include 'mach.p'
      integer i
      include 'caver.cm'
c this subroutine initializes and zeroes averages
      nblock=0
      do 1 i=1,maver
      avsum(i)=0.0
      avsq(i)=0.0
      ansum(i)=0.0
      ablock(i)=0.0
1     continue
      return
      end
