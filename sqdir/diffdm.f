      program diffdm
      implicit real*8 (a-h,o-z)
      integer mnr,mx,mnprm,ipickoff,intread,ip,ix,ist,nxs,isize
      integer nr,ifile,ndata,i,nfiles,mmin,nskip
      integer ln,i1,koff,l,koffi,k,n,iplot,igrid,ifgr,ind
      integer idim
      real*8 x,y,xv,rv,a1,a2,a3,a4
      parameter (mnr=6,mx=10000,mnprm=12)
      dimension idim(mnr),x(mx),ind(mnr)
     +,ifgr(mnr),rv(mx,2),xv(mx,2),ist(2),nxs(2),isize(2)
     +,y(mx,2)
      character  p(mnprm)*28,filen*14,ft*14,grid(mnprm,2)*28
      logical ifex

c reads a standard file and creates special graphics files
      mmin=mx
      nfiles=0
1     nfiles=nfiles+1
      if(nfiles.eq.3) go to 2

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

20    open(1,file=filen(1:ln),status='old',form='formatted')
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
      elseif(p(1).eq.'GRID'.and.intread(p(2)).eq.1) then
c save grid
         isize(nfiles)=idim(1)
         if(idim(1).lt.mmin) then
                mmin=idim(1)
                igrid=nfiles
         endif
         do i=2,n
         grid(i-1,nfiles)=p(i)
         enddo
         call setgrid(mx,isize(nfiles),rv(1,nfiles),p(3))
      endif
      go to 50
60    write (*,*)'*************end of record***'
      if(ifile.le.nskip) go to 100
      close(1)

      iplot=1
      write (*,*)' input all starting indices '
      read (*,*) (ind(i),i=1,nr)
      ind(1)=1

c compute offsets; i1 is first piece of data, koffi is spacing 
       koff=1
       i1=1
       do i=1,nr
         i1=i1+(ind(i)-1)*koff
         if(i.eq.iplot) koffi=koff
         koff=koff*idim(i)
       enddo

c move data from x to xv
       k=1
       do l=ind(iplot),idim(iplot)
         xv(k,nfiles)=x(i1)
         k=k+1
         i1=i1+koffi
       enddo
       ist(nfiles)=ind(iplot)
       nxs(nfiles)=idim(iplot)-ind(iplot)+1
       go to 1

2      continue

c subtract them
c set  grid to larger one
       ip=3-igrid
         call setgrid(mx,isize(ip),rv(1,ip),grid(2,ip))
       do i=1,nxs(igrid)
         call interp(rv(i,igrid),ix,a1,a2,a3,a4)
      y(i,1)=a1*xv(ix-1,ip)+a2*xv(ix,ip)+a3*xv(ix+1,ip)+a4*xv(ix+2,ip)
      y(i,2)=xv(i,igrid)
      xv(i,1)=y(i,1)-y(i,2)
         write (90,90) rv(i,igrid),y(i,1),y(i,2),xv(i,1)
90       format(4e15.6)
       enddo

       stop
       end
