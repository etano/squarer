      function intread(p)
      include 'mach.p'
      integer intread,ln
      character p*(*),blank*16,both*16
      data blank/'                '/
      ln=index(p,' ')-1
      both=blank(1:16-ln)//p(1:ln)
      read(both,1) intread
1     format(i16)
      return
      end
