subroutine constrnt(neq,nhbw,ndf,ncon,icon,vcon,glk,glm,glf,trm,mxneq)
include 'common.h'

integer, intent(in)           :: neq,nhbw,ndf,ncon,mxneq, icon(9)
real (kind=dbl) ,intent(in)   :: vcon(9)
real (kind=dbl) ,intent(inout):: trm(mxneq,mxneq), glm(mxneq,mxneq),glk(mxneq,mxneq),glf(mxneq)

!           ____________________________________________________________________
!
!           the subroutine is called in main to implement specified constraint
!            conditions (e.g., inclined supports) on the condensed system of
!            equations. array glm is used here as a temporary storage array.
!           ____________________________________________________________________

  !pi=3.14159265d0
  real (kind=dbl), parameter :: pi = dacos(-1.d0)
  real (kind=dbl) :: beta
  integer         :: i,j,k,jc,l,ic, idof


! include specified constraint conditions

  do ic=1,neq
    do jc=1,neq
      glm(ic,jc)=0.0d0
      trm(ic,jc)=0.0d0
    enddo
    trm(ic,ic)=1.0d0
  enddo

  do ic=1,ncon
    beta=vcon(ic)*pi/180.0d0
    idof=ndf*icon(ic)-1
    trm(idof,idof)    = dcos(beta)
    trm(idof,idof+1) = dsin(beta)
    trm(idof+1,idof) =-dsin(beta)
    trm(idof+1,idof+1)= dcos(beta)
  enddo

  l=0
  do  i=1,neq
    do j=1,nhbw
      glm(i,l+j)=glk(i,j)
    enddo
    l=l+1
  enddo

  do i=1,neq
    do j=i,neq
      glm(j,i)=glm(i,j)
    enddo
  enddo

  do i=1,neq
    do j=1,neq
      glk(i,j)=glm(i,j)
    enddo
  enddo

  do i=1,neq
    do j=1,neq
      glm(i,j)=0.0d0
      do k=1,neq
        glm(i,j)=glm(i,j)+trm(i,k)*glk(k,j)
      enddo
    enddo
  enddo

  do i=1,neq
    do j=1,neq
      glk(i,j)=0.0d0
      do k=1,neq
        glk(i,j)=glk(i,j)+glm(i,k)*trm(j,k)
      enddo
    enddo
  enddo

  do i=1,neq
    do j=1,neq
      trm(i,j)=glk(i,j)
    enddo
  enddo

  l=0
  do i=1,neq
    do j=1,nhbw
      glk(i,j)=trm(i,l+j)
    enddo
    l=l+1
  enddo

  do i=1,neq
    glm(i,1)=0.0d0
    do k=1,neq
      glm(i,1)=glm(i,1)+trm(i,k)*glf(k)
    enddo
      glf(i)=glm(i,1)
  enddo

  return
end ! end of subroutine constrnt() !
