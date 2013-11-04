c        program main
c        call test
c        end
 
        subroutine test
        implicit none
        integer m,maxm,i,j,jx
        parameter (maxm=6)
        double precision s(0:maxm,0:2*maxm+1)
        double precision p(0:2*maxm+1)
        double precision x,dx
        write(*,*) 'Insert m'
        read(*,*) m
        call basis(m,maxm,s)
        do i=0,m
         write(unit=*,fmt=50) (s(i,j),j=0,2*m+1)
        enddo
        dx=0.01
        do jx=0,100
         x=dx*dble(jx)
         do i=0,m
          p(i)=0.d0
          do j=0,2*m+1
           p(i)=p(i)+s(i,j)*x**j
          enddo
         enddo
         write(unit=20,fmt=50) x,(p(i),i=0,m)
        enddo
50      format(20(f14.6,1x))
        end

        subroutine basis(m,maxm,s)
        implicit none
        integer m,maxm,nd
        integer ialpha,ien,i,j,ibeta,k
        integer info,job 
        integer fact(0:2*m+1)
        double precision work(m+1),ipiv(m+1),det(2)
        double precision a(m+1,m+1),b(m+1,m+1)
        double precision x(m+1,m+1)
        integer delta(0:m,0:m)
        double precision s(0:maxm,0:2*maxm+1)
 
        nd=m+1

!    \\ LOAD FACTORIAL
        fact(0)=1
        do i=1,2*m+1
         fact(i)=fact(i-1)*i
        enddo 

!    \\ LOAD MATRIX A
        do ien=1,nd
         do ialpha=1,nd
          a(ien,ialpha)=dble(fact(ien+m))/dble(fact(ien+m-ialpha+1))
         enddo
        enddo

!    \\ LOAD MATRIX B
        do ialpha=1,nd
         do ibeta=1,nd
          if (ialpha.gt.ibeta) then
           b(ibeta,ialpha)=0.d0
          else
           b(ibeta,ialpha)=-1.d0/dble(fact(ibeta-ialpha))
          endif
         enddo
        enddo

!    \\ INVERT MATRIX A
!       call dgetrf(nd,nd,a,nd,ipiv,info)
!       call dgetri(nd,a,nd,ipiv,work,nd,info)
        call sgefa(a,nd,nd,ipiv,info)
        job=1
        call sgedi(a,nd,nd,ipiv,det,work,job)

!    \\ COMPUTE MATRIX X
        do i=1,nd
         do j=1,nd
          x(i,j)=0.d0
          do k=1,nd
           x(i,j)=x(i,j)+b(i,k)*a(k,j)
          enddo
         enddo
        enddo

!    // LOAD UNIT MATRIX
        do i=0,m
         do j=0,m
          delta(i,j)=0
         enddo
         delta(i,i)=1
        enddo

!    // LOAD FINAL MATRIX
        do i=0,m
         do j=0,m
          s(i,j)=dble(delta(i,j))/dble(fact(i))
          s(i,j+m+1)=x(i+1,j+1)
         enddo
        enddo

        return
        end 
