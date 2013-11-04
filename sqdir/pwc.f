      subroutine pwc(ndim,ldim,cl)
c returns the partial wave normalization in 1,2,3 dimensions
      real*8 cl(0:ldim),pii
        pii=.31830 98861 83790d0
        if(ndim.eq.1) then
         cl(0)=1.
         if(ldim.gt.0) then
             write(*,*) "ldim=",ldim
!            call warning("ldim > 0")
         endif
        elseif(ndim.eq.2) then
         cl(0)=0.5d0*pii
         do l=1,ldim
           cl(l)=pii
         enddo
        elseif(ndim.eq.3) then
         do l=0,ldim
           cl(l)=0.25d0*(2*l+1)*pii
         enddo
        else
          write (*,*) 'error in pwc. ndim=',ndim
          stop
        endif
        return
        end
