      function convert(type,from,to)
      implicit none
c convert units-only energy and length work
c to add units. put in the name, its value in mks and bump nu by 1.
      integer nunits,nspaces,n,i,j
      parameter (nunits=2,nspaces=10) 
      character type*(*),from*(*),to*(*),name(nspaces,nunits)*4
     +,p*4 
      integer nu(nunits) 
      real*8  size(nspaces,nunits) ,x(2) , convert
      save nu,size,name 
c energy
      data nu(1)/5/
      data name(1,1)/'J'/
      data size(1,1)/1./
      data name(2,1)/'eV'/
      data size(2,1)/1.60217733e-19/
      data name(3,1)/'R'/
      data size(3,1)/2.1798741e-18/
      data name(4,1)/'K'/
      data size(4,1)/1.380658e-23/
      data name(5,1)/'H'/
      data size(5,1)/4.3597482e-18/
c length
      data nu(2)/4/
      data name(1,2)/'m'/
      data size(1,2)/1./
      data name(2,2)/'nm'/
      data size(2,2)/1.e-9/
      data name(3,2)/'A'/
      data size(3,2)/1.e-10/
      data name(4,2)/'a0'/
      data size(4,2)/.529177249e-10/
 
      if(type.eq.'energy') then
        n=1
      elseif(type.eq.'length') then
        n=2
      else
        write (*,*)'Convert cant do ',type
        stop
      endif
 
c find ratios in mks
       do j=1,2    
       x(j)=0. 
       if(j.eq.1)p=from  
       if(j.eq.2)p=to 
 
       do i=1,nu(n)
       if(name(i,n).eq.p) x(j)=size(i,n)
       enddo
 
       if(x(j).le.0) then
          write (*,*)'Convert doesnt know ',p
          stop
       endif
 
       enddo

       convert=x(1)/x(2)
       return
       end
