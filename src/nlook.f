      program nlook
      implicit real*8 (a-h,o-z)
      parameter (mnr=6,mx=20000,mnprm=12,ngm=1000)
      dimension idim(mnr),x(mx),ind(mnr)
     +,ifgr(mnr),rv(ngm,mnr),xv(ngm,20)
      character  p(mnprm)*28,filen*14,ft*14
      logical ifex

c reads a standard file and creates special graphics files
10    write (*,*) ' what file? stop to stop'
      read (*,9) ft
9     format(a14)
      if(ft.eq.'!!')then
        ft=filen
      elseif(ft.eq.'stop') then
        stop
      else
        filen=ft
      endif

      ln=index(filen,' ')-1
c does this file exist?
      inquire(file=filen(1:ln),exist=ifex)
      if(ifex) go to 20
      write (*,*)' file does not exist '
      go to 10

20    open(1,file=filen(1:ln),status='old',form='formatted',err=8010)
      write (*,*) ' how many files to skip ?'
      read (*,*) nskip
      ifile=0
100   ifile=ifile+1
       nr=0
      do 30 i=1,mnr
30    idim(i)=0
c get rank nr and the number of components for each rank tensor-index
      call checkrnk(1,nr,idim)
      ndata=1
      do 40 i=1,nr
      ifgr(i)=0
40    ndata=ndata*idim(i)
      if(ndata.gt.mx) then
        write (*,*)' dimension of x too small '
        stop
      endif

50    if(ipickoff(1,p,n,mnprm).ne.0.or.n.le.0) then
        write (*,*)' data file empty ',n
        stop
      endif

      if(p(1).eq.'BEGIN') then
        read (1,*) (x(k),k=1,ndata)
        go to 60 
      elseif(p(1).eq.'GRID') then
        j=intread(p(2))
        ifgr(j)=1
        if(p(3).eq.'DISCRETE') then
c read the discrete grid points
          read(1,*,err=8020) (rv(i,j),i=1,idim(j))
        else
c the grid is either LINEAR or LOG
          call setgrid(mx,idim(j),rv(1,j),p(3))
        endif
      endif
      go to 50
60    write (*,*)'*************end of record***'
      if(ifile.le.nskip) go to 100

c we plot along the first axis
      iplot=1
      nplots=0

70    write (*,*)' input all starting indices (a 0 to stop)'
      read (*,*) (ind(i),i=1,nr)
      if(ind(1).le.0) go to 77
      nplots=nplots+1
      ind(iplot)=1
      
c compute offsets; i1 is first piece of data, koffi is spacing 
       koff=1
       i1=1
       do 80 i=1,nr
       i1=i1+(ind(i)-1)*koff
       if(i.eq.iplot) koffi=koff
80     koff=koff*idim(i)

c move data from x to xv
       k=1
       do 90 l=ind(iplot),idim(iplot)
       xv(k,nplots)=x(i1)
       k=k+1
90     i1=i1+koffi
       i1=ind(iplot)
       nxs=idim(iplot)-i1+1
       go to 70
c finshed
77     continue
      ind(iplot)=1
        write (*,*) nxs,nplots
       do i=1,nxs
       write (99,99) rv(ind(iplot)+i-1,iplot),(xv(i,l),l=1,nplots)
99     format(12e15.7)
       enddo
       stop

8010   write(*,*) ' flook: can not open file ',filen(1:ln)
       stop
8020   write(*,*) ' flook: can not read file ',filen(1:ln)
       stop
       end
