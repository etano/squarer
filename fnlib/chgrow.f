       subroutine chgrow(ratinv,nrow,rewrow,irow,rowi,dot,di)
      include 'mach.p'
cc       update a row of the inverse matrix
c   nrow is order,   irow is column number
c    ratinv is old inverse,   rewrow is replacement row
c this is the Sherman-Morrison update formula
       real*8 ratinv,rewrow,rowi,dot,di,f,t
       integer nrow,irow,i,k,j
       dimension ratinv(nrow,nrow),rewrow(nrow),rowi(nrow),dot(nrow)
       if(nrow.le.0) return
       f=-1.d0/di
       do 1 i=1,nrow
       dot(i)=0.d0
 1     rowi(i)=ratinv(irow,i)
       do 12 k=1,nrow
       t=f*rewrow(k)
       do 12 i=1,nrow
 12    dot(i)=dot(i)+t*ratinv(i,k)
       do 3 j=1,nrow
       do 3 i=1,nrow
3      ratinv(i,j)=ratinv(i,j)+dot(i)*rowi(j)
       do 20 i=1,nrow
20     ratinv(irow,i)=-f*rowi(i)
       return
       end
