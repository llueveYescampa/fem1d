subroutine eqnsolvr(nrm,ncm,neqns,nbw,band,rhs,ires)
include 'common.h'

integer, intent(in)           :: nrm,ncm,neqns,nbw,ires
real (kind=dbl) ,intent(inout):: band(nrm,ncm), rhs(nrm)

! _________________________________________________________________

! the subroutine is called in main to solve symmetric and banded set
! of equations using the gauss elimination method:[band]{u} = {rhs}.
! the coefficient matrix is input as band(neqns,nbw) and the column
! vector is input as rhs(neqns), where neqns is the actual number
! of equations and nbw is the half band width. the true dimensions
! of the matrix [band] in the calling program, are nrm by ncm. when
! ires is greater than zero, the right hand elimination is skipped.
! _________________________________________________________________


  real (kind=dbl) :: factor
  integer         :: icol, ijk, jcol, jki, lstsub, npiv, npivot, nrow, meqns, ncol

  meqns=neqns-1
  if(ires.le.0) then
    do npiv=1,meqns
      npivot=npiv+1
      lstsub=npiv+nbw-1
      if(lstsub.gt.neqns) then
        lstsub=neqns
      endif

      do nrow=npivot,lstsub
        ncol=nrow-npiv+1
        factor=band(npiv,ncol)/band(npiv,1)
        do ncol=nrow,lstsub
          icol=ncol-nrow+1
          jcol=ncol-npiv+1
          band(nrow,icol)=band(nrow,icol)-factor*band(npiv,jcol)
        enddo
        rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
      enddo
    enddo

  else
    do npiv=1,meqns
      npivot=npiv+1
      lstsub=npiv+nbw-1
      if(lstsub.gt.neqns) then
        lstsub=neqns
      endif
      do nrow=npivot,lstsub
        ncol=nrow-npiv+1
        factor=band(npiv,ncol)/band(npiv,1)
        rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
      enddo
    enddo
  endif

! back substitution

  do ijk=2,neqns
    npiv=neqns-ijk+2
    rhs(npiv)=rhs(npiv)/band(npiv,1)
    lstsub=npiv-nbw+1
    if(lstsub.lt.1) then
      lstsub=1
    endif
    npivot=npiv-1
    do jki=lstsub,npivot
      nrow=npivot-jki+lstsub
      ncol=npiv-nrow+1
      factor=band(nrow,ncol)
      rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
    enddo
  enddo
  rhs(1)=rhs(1)/band(1,1)
  return
end ! end of subroutine eqnsolvr() !
