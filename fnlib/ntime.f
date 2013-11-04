      function ntime(i,iref,irev,nslices)
      include 'mach.p'
      integer ntime,i,iref,irev,nslices
      if(iref.ne.0) then
         ntime=i-iref
         if(ntime.le.0)ntime=ntime+nslices
         if(ntime.gt.irev)ntime=ntime-nslices
      else
         ntime=nslices
      endif
      return
      end
