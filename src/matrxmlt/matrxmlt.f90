subroutine matrxmlt(mxneq,n,a,b,c)
include 'common.h'

integer, intent(in)          :: mxneq,n
real (kind=dbl) ,intent(in)  :: a(mxneq,mxneq),b(mxneq,mxneq)
real (kind=dbl) ,intent(out) :: c(mxneq,mxneq)

! ________________________________________________________________
! Called in EGNSOLVR to compute the product of matrices [A] & [B]:
!                         [C]=[A][B]
! ________________________________________________________________

  integer :: i,j,k
  do i=1,n
    do j=1,n
      c(i,j)=0.0
      do k=1,n
        c(i,j)=c(i,j)+a(i,k)*b(k,j)
      enddo
    enddo
  enddo

  return
end ! end of subroutine matrxmlt() !
