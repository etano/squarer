       subroutine cutkset(ndim,vol,a,alpha,itol,tau,cut)
       include 'mach.p'
       integer itol,iter,ndim,icon,i
       real*8 a(ndim),all,alpha,tol,ak0,ak1,f0,fk,cut,tau,vol
       parameter (iter=1000)
   
       icon = 0
       all = 0.25/alpha**2
       tol = 10.**(-itol)*vol/12./tau
       ak0 = a(1)
       do i=1,ndim
        ak0 = min(ak0,a(i))
       enddo
       f0 = exp(-ak0*ak0*all)/ak0/ak0
       ak1 = ak0
1      icon = icon + 1
       if (icon.gt.iter) then
        write (*,*) 'subroutine setcutk : more than ',iter,' iterations'
        stop
       end if
       ak1 = ak1 + ak0
       fk = exp(-ak1*ak1*all)/ak1/ak1
       if (fk.gt.tol) goto 1
       cut = ak1 - ak0
       fk = exp(-cut*cut*all)/cut/cut
       write (*,*)
       write (*,*)
       write (*,*) 'subroutine setcutk'
       write (*,*)
       write (*,*) '|k|(max),  fk/f0 : ',cut,fk
       return
       end
