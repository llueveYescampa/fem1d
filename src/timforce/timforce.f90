subroutine timforce(elf,elx,fx0,fx1,fx2,h,ntype,ne,f3,mxelm)
include 'common.h'

integer, intent(in)           :: ntype,ne,mxelm
real (kind=dbl) ,intent(in)   :: elx(4), fx0,fx1,fx2,h
real (kind=dbl) ,intent(inout)  :: elf(9), f3(mxelm)

!        __________________________________________________________________

!           called in coeffcnt to compute element force vector for the
!                consistent interpolation timoshenko element (cie)
!        __________________________________________________________________

  include 'shp.h'
!     dimension gauspt(5,5),gauswt(5,5),elf(9),elx(4),ex(3),f3(mxelm)



  real (kind=dbl) :: xi, const, fx, ex(3), x, gauspt(5,5),gauswt(5,5)
  integer         :: i, j, npe,ni, ndf, ii, iel, ngp

  data gauspt/5*0.0d0,-.57735027d0,.57735027d0,3*0.0d0,-.77459667d0, &
    0.0d0,.77459667d0,2*0.0d0,-.86113631d0,-.33998104d0,.33998104d0,  &
    .86113631d0,0.0d0,-.906180d0,-.538469d0,0.0d0,.538469d0,.906180d0/

  data gauswt/2.0d0,4*0.0d0,2*1.0d0,3*0.0d0,.55555555d0,.88888888d0, &
    0.55555555d0,2*0.0d0,.34785485d0,2*.65214515d0,.34785485d0,0.0d0, &
    0.236927d0,.478629d0,.568889d0,.478629d0,.236927d0/


  npe=3
  iel=2
  ndf=2
  ngp=iel+1
  do  i=1,6
    elf(i)=0.0d0
  end do

  ex(1)=elx(1)
  ex(2)=elx(1)+0.5*h
  ex(3)=elx(2)

  do ni=1,ngp
    xi=gauspt(ni,ngp)
    call shape1d(h,iel,npe,xi)
    const=gj*gauswt(ni,ngp)
    x = 0.0d0
    do  j=1,npe
      x = x + sf(j)*ex(j)
    end do

!   compute the polynomial variation of fx

    if(ntype.eq.2)then
      fx=fx0+(fx1+fx2*x)*x
    else
      fx=(fx0+fx1*x)*x
    endif

!   element force vector for the consistent interpolation beam element

    ii=1
    do  i=1,npe
      elf(ii)=elf(ii)+fx*sf(i)*const
      ii=ndf*i+1
    end do
  end do

! rearrange the element coefficients

  f3(ne)=elf(3)
  elf(1)=elf(1)+0.5*f3(ne)
  elf(2)=-0.125*f3(ne)*h
  elf(3)=elf(5)+0.5*f3(ne)
  elf(4)= 0.125*f3(ne)*h
  return
end ! end of subroutine timforce() !
