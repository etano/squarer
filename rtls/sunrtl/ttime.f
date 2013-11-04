      integer i(3),n,time
          call itime(i)
          write (*,*) ' hr:min:sec',i
      call idate(i)
          write (*,*) ' day mn yr ',i
      n=time()
      write (*,*) ' time since 1970 ',n
      stop
        end
