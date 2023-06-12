subroutine assemble(nod,mxelm,mxneq,ndf,npe,ne,item,glk,glm,glf)
include 'common.h'

integer, intent(in)            :: nod(mxelm,4), mxelm,mxneq,ndf,npe,ne,item
real (kind=dbl) ,intent(inout) :: glf(mxneq), glk(mxneq,mxneq), glm(mxneq,mxneq)

!       __________________________________________________________________

!         the subroutine is called in main to assemble element coefficient
!         matrices (in a upper banded matrix form) and right-hand vectors

!          {elf}.... element source vector, {f}
!          {elk}.... element coefficient matrix, [k]
!          {elm}.... element coefficient matrix, [m]
!          [nod].... connectivity matrix, [b]
!       __________________________________________________________________

  include 'stf1.h'

  integer   :: i, ii, j, jj, l, m, nc, ncl, nr


  if(item <= 2) then

!   assemble element coefficient matrix elk and source vector elf

    do i = 1, npe
      nr = (nod(ne,i) - 1)*ndf
      do ii = 1, ndf
        nr = nr + 1
        l = (i-1)*ndf + ii
        glf(nr) = glf(nr) + elf(l)
        do j = 1, npe
          ncl = (nod(ne,j)-1)*ndf
          do jj = 1, ndf
            m = (j-1)*ndf + jj
            nc = ncl-nr+jj+1
            if(nc > 0 ) then ! 20,20,10
              glk(nr,nc) = glk(nr,nc) + elk(l,m)
            endif
          enddo
        enddo
      enddo
    enddo
  else

! assemble element matrices into full global matrices

    do i=1,npe
      nr=(nod(ne,i)-1)*ndf
      do ii=1,ndf
        nr=nr+1
        l=(i-1)*ndf+ii
        do j=1,npe
          nc=(nod(ne,j)-1)*ndf
          do jj=1,ndf
            m=(j-1)*ndf+jj
            nc=nc+1
            glk(nr,nc)=glk(nr,nc)+elk(l,m)
            glm(nr,nc)=glm(nr,nc)+elm(l,m)
          enddo
        enddo
      enddo
    enddo

  endif
  return
end ! end of subroutine assemble() !
