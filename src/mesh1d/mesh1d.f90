subroutine mesh1d(nem,npe,nod,mxelm,mxnod,dx,glx)
include 'common.h'

integer, intent(in)           :: nem,npe,mxelm,mxnod
integer, intent(out)          :: nod(mxelm,4)
real (kind=dbl) ,intent(in)   :: dx(mxnod)
real (kind=dbl) ,intent(inout):: glx(mxnod)
!       __________________________________________________________________

!       The subroutine is called in MAIN to compute arrays {GLX} and [NOD]

!          {GLX}.... Vector of global coordinates
!          {DX}..... Vector of element lengths [DX(1) = node 1 coordinate]
!          [NOD].... Connectivity matrix
!       __________________________________________________________________

  integer :: i, ii, n

  !       Generate the elements of the connectivity matrix
  do i=1,npe
    nod(1,i)=i
  end do
  do n=2,nem
    do i=1,npe
      nod(n,i) = nod(n-1,i)+npe-1
    end do
  end do

  !       Generate global coordinates of the global nodes

  glx(1)=dx(1)
  if(npe.eq.2)then
    do  i=1,nem
      glx(i+1) = glx(i) + dx(i+1)
    end do
  else
    do i=1,nem
    ii=2*i
      glx(ii) = glx(ii-1) + 0.5*dx(i+1)
      glx(ii+1)=glx(ii-1) + dx(i+1)
    end do
  endif
  return
end ! end of subroutine mesh1d() !
