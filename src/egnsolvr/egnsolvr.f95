subroutine egnsolvr(n,a,b,xx,x,negn,nr,mxneq)
include 'common.h'

integer, intent(in)             :: n,negn,mxneq
integer, intent(out)            :: nr
real (kind=dbl) ,intent(in)     :: a(mxneq,mxneq),x(mxneq,mxneq)
real (kind=dbl) ,intent(inout)  :: b(mxneq,mxneq),xx(mxneq)

!__________________________________________________________________
!
! the subroutine is called in main to solve the eigenvalue problem
!
!                     [a]{x} = lambda[b]{x}
!
! the program can be used only for positive-definite [b] matrix.
! the dimensions of v, vt, w, and ih should be equal to mxneq.
!__________________________________________________________________


  real (kind=dbl), dimension(:,:), allocatable :: v
  real (kind=dbl), dimension(:,:), allocatable :: vt
  real (kind=dbl), dimension(:,:), allocatable :: w
  integer        , dimension(:)  , allocatable :: ih

  !real (kind=dbl) :: v(500,500),vt(500,500),w(500,500)
  integer         :: i, j ! , ih(500)


  allocate ( v(mxneq,mxneq))
  allocate ( vt(mxneq,mxneq))
  allocate ( w(mxneq,mxneq))
  allocate ( ih(mxneq) )

! call subroutine jacobi to diagonalize [b]

  call jacobi (n,b,negn,nr,v,xx,ih,mxneq)

! make diagonalized [b] symmetric

  do  i=1,n
    do  j=1,n
       b(j,i)=b(i,j)
    enddo
  enddo

! check (to make sure) that [b] is positive-definite

  do i=1,n
    if (  b(i,i) < 0.0 ) then !  20,30,30
      write(6,"(/'*** matrix [glm] is not positive-definite ***')")
      stop
    endif
  enddo

!       the eigenvectors of [b] are stored in array v(i,j)
!       form the transpose of [v] as [vt]

  do i=1,n
    do j=1,n
      vt(i,j)=v(j,i)
    enddo
  enddo


! find the product [f]=[vt][a][v] and store in [a] to save storage

  call matrxmlt (mxneq,n,vt,a,w)
  call matrxmlt (mxneq,n,w,v,a)

! get [gi] from diagonalized [b], but store it in [b]

  do i=1,n
    b(i,i)=1.0/dsqrt(b(i,i))
  enddo

! find the product [q]=[gi][f][gi]=[b][a][b] and store in [a]

  call matrxmlt (mxneq,n,b,a,w)
  call matrxmlt (mxneq,n,w,b,a)

! we now have the form [q]{z}=lamda{z}. diagonalize [q] to obtain
! the eigenvalues by calling jacobi. the eigenvalues are returned
! as diag [a].

  call jacobi (n,a,negn,nr,vt,xx,ih,mxneq)
  do j=1,n
    xx(j)=a(j,j)
  enddo

!       the eigenvectors are computed from the relation,
!                       {x}=[v][gi]{z}=[v][b][vt]
!       since {z} is stored in [vt]

  call matrxmlt (mxneq,n,v,b,w)
  call matrxmlt (mxneq,n,w,vt,x)


  deallocate ( v )
  deallocate ( vt )
  deallocate ( w )
  deallocate ( ih )

  return
end ! end of subroutine egnsolvr() !
