subroutine echodata(in,it)
Use, intrinsic :: iso_fortran_env, Only : iostat_end
include 'common.h'

integer, intent(in) :: in, it

  character (len=80) :: aa
  integer :: ioerr

  write(it,"(5x,'*** echo of the input data starts ***'/)")

  do  !
    read(in,'(80a)', iostat = ioerr) aa
    if ( ioerr == iostat_end  ) exit
    write(it,'(1x,80a)') aa
  enddo ! while

  write(it,"(5x,'**** echo of the input data ends ****'/)")
  rewind(in)
  return
end ! end of subroutine echodata() !
