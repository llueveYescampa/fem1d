subroutine postproc(dcax,dcbx,dccx,f3,glf,glx,nod,icont,ielem,npe,model,ntype,item,mxelm,mxneq,mxnod,nem,ndf)
include 'common.h'

integer, intent(in)            :: nod(mxelm,4),icont,ielem,npe,model,ntype,item,mxelm,mxneq,mxnod,nem,ndf
real (kind=dbl) ,intent(inout) :: dcax(mxelm,2),dcbx(mxelm,2),dccx(mxelm,2), f3(mxelm),glf(mxneq),glx(mxnod)

!       __________________________________________________________________

!       the subroutine is called in main to compute the solution and its
!       derivatives at five points, including the nodes of the element.
!       the bending moment (bm) and shear force (vf) are computed as per
!       the definitions given in fig. 5.2.1 and eq. (5.3.1) of the book;

!           x........ global (i.e., problem) coordinate
!           xi ...... local (i.e., element) coordinate
!           sf....... element interpolation (or shape) functions
!           gdsf..... first derivative of sf w.r.t. global coordinate
!           gddsf.... second derivative of sf w.r.t. global coordinate
!           elu...... element solution vector
!           u........ interpolated solution
!           du....... interpolated derivative of the solution
!           w........ interpolated transverse deflection
!           s........ interpolated rotation function
!           ds....... interpolated derivative of the rotation
!           dw....... interpolated derivative of the transverse deflection
!           ddw...... interpolated second derivative of trans. deflection
!           dddw..... interpolated third derivative of trans. deflection
!        __________________________________________________________________

  real (kind=dbl) :: xp(9),elx(4),elu(9)
  real (kind=dbl) :: anu21,u,vf,w,xc,xi,sr,st,theta,sv,sfv,psi,bm,bmr,bmt,c11,c12,c22,d11,d12,d22,ddw
  real (kind=dbl) :: dddw,denom,di,dsi,du,dw,h,dpsi=0.0d0,si
  integer         :: npts,i,j,l,li,ne,ni

  include 'io.h'
  include 'shp.h'
  include 'stf2.h'

  data xp/-1.0d0, -0.750d0, -0.50d0, -0.250d0, 0.0d0, 0.250d0,0.50d0, 0.750d0, 1.0d0/

  npts=9
  do ne = 1, nem
    if(icont /= 1) then
      ax0=dcax(ne,1)
      ax1=dcax(ne,2)
      bx0=dcbx(ne,1)
      bx1=dcbx(ne,2)
      cx0=dccx(ne,1)
      cx1=dccx(ne,2)
    endif
    l=0
    do i=1,npe
      ni=nod(ne,i)
      if(icont /= 1)then
        elx(1)=0.0
        elx(2)=0.5*glx(ne)
        elx(npe)=glx(ne)
      else
        elx(i)=glx(ni)
      endif
      li=(ni-1)*ndf
      do j=1,ndf
        li=li+1
        l=l+1
        elu(l)=glf(li)
      enddo
    enddo
    h = elx(npe) - elx(1)
    do ni=1,npts
      xi = xp(ni)
      call shape1d(h,ielem,npe,xi)
      if(model == 3)then
        w=0.0
        dw=0.0
        ddw=0.0
        xc=elx(1)+0.5*h*(1.0+xi)
        do i=1,4
          w =w + sf(i)*elu(i)
          dw =dw + gdsf(i)*elu(i)
          ddw=ddw+ gddsf(i)*elu(i)
        enddo
        dddw=((elu(1)-elu(3))*2.0/h-(elu(4)+elu(2)))*6.0/(h*h)
        theta=-dw
        if(ntype == 0)then
          bm=-(bx0+xc*bx1)*ddw
          vf=-(bx0+xc*bx1)*dddw - bx1*ddw
          write(it,'(2x,6e13.5)')xc,w,theta,bm,vf
        else
          anu21=bx0*ax0/ax1
          di=(bx1**3)/12.0
          d11=di*ax0/(1.0-bx0*anu21)
          d22=d11*(ax1/ax0)
          d12=bx0*d22
          bmr=-(d11*ddw*xc+d12*dw)
          bmt=-(d12*ddw*xc+d22*dw)
          if(xc  /=  0.0d0)then
            sfv=-d11*(xc*dddw+ddw)+d22*dw/xc
            write(it,'(2x,6e13.5)')xc,w,theta,bmr,bmt,sfv
          else
            write(it,'(2x,6e13.5)')xc,w,theta,bmr,bmt
          endif
        endif
      else
        xc=0.0
        do i=1,npe
          xc=xc+sf(i)*elx(i)
        enddo
        if(model == 1)then
        u=0.0
        du=0.0
        do i=1,npe
        u=u+sf(i)*elu(i)
        du=du+gdsf(i)*elu(i)
        enddo
        if(ntype == 0)then
          sv=(ax0+ax1*xc)*du
          write(it,'(2x,6e13.5)')xc,u,sv
        else
          anu21=bx0*ax0/ax1
          if(ntype == 1)then
            c11=bx1*ax0/(1.0-bx0*anu21)
            c22=c11*(ax1/ax0)
            c12=bx0*c22
          else
            denom=1.0-bx0-anu21
            c11=bx1*ax0*(1.0-bx0)/(1.0+bx0)/denom
            c22=bx1*ax1*(1.0-anu21)/(1.0+anu21)/denom
            c12=bx0*c22
          endif
          if(xc /= 0.0d0)then
            sr=c11*du+c12*u/xc
            st=c12*du+c22*u/xc
            write(it,'(2x,6e13.5)')xc,u,sr,st
          else
            write(it,'(2x,6e13.5)')xc,u,du
          endif
        endif
        else
          !       model == 2   calculations
          if(item == 0 .and. ntype.gt.1)then
            h=elx(npe)-elx(1)
            call timstres(ax0,elu,xi,w,dw,si,dsi,ne,f3,h,mxelm)
          else
            w   =0.0
            dw =0.0
            psi =0.0
            dpsi =0.0
            do i=1,npe
              l=2*i-1
              w   = w   +sf(i)*elu(l)
              dw = dw +gdsf(i)*elu(l)
              psi = psi +sf(i)*elu(l+1)
              dpsi= dpsi+gdsf(i)*elu(l+1)
            enddo
          endif
          if(ntype == 0 .or. ntype == 2)then
            bm=(bx0+bx1*xc)*dpsi
            vf=(ax0+ax1*xc)*(dw+psi)
            write(it,'(2x,6e13.5)')xc,w,psi,bm,vf
          else
            anu21=bx0*ax0/ax1
            di =(bx1**3)/12.0
            d11=di*ax0/(1.0-bx0*anu21)
            d22=d11*(ax1/ax0)
            d12=bx0*d22
            bmr=(d11*dpsi*xc+d12*psi)
            bmt=(d12*dpsi*xc+d22*psi)
            sfv=fx2*(dw+psi)*xc
            write(it,'(2x,6e13.5)')xc,w,psi,bmr,bmt,sfv
          endif
        endif
      endif
    enddo
  enddo
  return
end ! end of subroutine postproc() !
