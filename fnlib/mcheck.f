      subroutine mcheck(iv,mv,civ,cmv,sub)
      include 'mach.p'
      integer iv,mv
      character civ*(*),cmv*(*),sub*(*)
      write (77,1) iv,mv,civ,cmv,sub
1     format(2i8,3a10)
      if(iv.le.mv) return
      write(*,2) sub,civ,cmv,iv,mv
2     format(' memory overflow in ',a,' variable ',a,' ',a,
     &' value ',i9,' limit ',i9)
      stop
      end
