subroutine timstres(ga,elu,xi,w,dw,s,ds,ne,f3,h,mxelm)
implicit none
integer, parameter :: dbl = selected_real_kind(p=13,r=200)

integer, intent(in)           :: ne,mxelm
real (kind=dbl) ,intent(in)   :: ga, xi,h, elu(9),f3(mxelm)
real (kind=dbl) ,intent(inout)  :: ds, s, w, dw


!     __________________________________________________________________

!     Called in POSTPROC to compute solution and its global derivatives
!       at nine points (including the nodes) of the Timoshenko element

!         XC........     Global (i.e., problem) coordinate
!         XI .......     Local (i.e., element) coordinate
!         SFL, SFQ..     Lagrange linear and quadratic shape functions
!         DSFL,DSFQ:     First derivative of SF w.r.t. global coordinate
!         ELU....... Column vector of generalized displacements
!         W, DW..... Transverse deflection and its derivative
!         S, DS..... Rotation and its derivative
!           __________________________________________________________________

!  implicit real*8 (a-h,o-z)

  include 'io.h'

!  dimension elu(9),sfl(2),sfq(3),dsfl(2),dsfq(3),f3(mxelm)

  real (kind=dbl) :: gj, w3
  real (kind=dbl) :: sfl(2),dsfl(2),dsfq(3), sfq(3)

  gj =       h*0.5

! Interpolation functions for the Lagrange LINEAR element

  sfl(1) = 0.5*(1.0-xi)
  sfl(2) = 0.5*(1.0+xi)
  dsfl(1) = -0.5/gj
  dsfl(2) = 0.5/gj

! Interpolation functions for the Lagrange QUADRATIC element

  sfq(1) = -0.5*xi*(1.0-xi)
  sfq(2) = 1.0-xi*xi
  sfq(3) = 0.5*xi*(1.0+xi)
  dsfq(1) = -0.5*(1.0-2.0*xi)/gj
  dsfq(2) = -2.0*xi/gj
  dsfq(3) = 0.5*(1.0+2.0*xi)/gj

  w3=(3.0*h*f3(ne)/ga + 8.0*(elu(1)+elu(3)) + 2.0*(elu(4)-elu(2))*h)/16.0
  w = sfq(1)*elu(1) + sfq(2)*w3 + sfq(3)*elu(3)
  dw= dsfq(1)*elu(1) +dsfq(2)*w3 +dsfq(3)*elu(3)
  s = sfl(1)*elu(2) + sfl(2)*elu(4)
  ds= dsfl(1)*elu(2) +dsfl(2)*elu(4)

  return
end ! end of subroutine timstres() !
