      subroutine chgclm(ratinv,ncolm,rewclm,icolm,rowi,dot,di)
      include 'mach.p'
c   update a column of the inverse matrix
c   ncolm is order,   icolm is column number
c    ratinv is old inverse,   rewclm is replacement column
c the Sherman-Morrison update formula
c
       real*8 ratinv,rewclm,rowi,dot,di,f,t
       integer ncolm,icolm,i,k,j
       dimension ratinv(ncolm,ncolm),rewclm(ncolm),rowi(ncolm)
     +  ,dot(ncolm)
       if(ncolm.le.0) return
       f=-1.d0/di
       do 1 i=1,ncolm
       dot(i)=0.d0
1      rowi(i)=ratinv(i,icolm)
       do 12 k=1,ncolm
       t=f*rewclm(k)
       do 12 i=1,ncolm
 12    dot(i)=dot(i)+t*ratinv(k,i)
       do 3 j=1,ncolm
       do 3 i=1,ncolm
3      ratinv(i,j)=ratinv(i,j)+dot(j)*rowi(i)
       do 20 i=1,ncolm
20     ratinv(i,icolm)=-f*rowi(i)
       return
       end
