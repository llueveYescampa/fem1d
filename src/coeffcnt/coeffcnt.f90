subroutine coeffcnt(ielem,item,model,ndf,npe,time,ntype,ne,f3,mxelm)
include 'common.h'

integer, intent(in)         :: ielem,item,model,ndf,npe,ntype,ne,mxelm
real (kind=dbl) ,intent(in) :: time,f3(mxelm)


!         __________________________________________________________________

!         the subroutine is called in main to compute coefficient matrices
!          and source vector for the model problem in eq. (1) (see main)

!             x......... global (i.e., problem) coordinate
!             xi ....... local (i.e., element) coordinate
!             h......... element length
!             {sf}...... element interpolation (or shape) functions
!             {gdsf}.... first derivative of sf w.r.t. x
!             {gddsf}... second derivative of sf w.r.t. x
!             gj........ determinant of the jacobian matrix
!             [gauspt].. 4x4 matrix of gauss points: n-th column corresponds
!                        to the n-point gauss rule
!            [gauswt].. 4x4 matrix of gauss weights (see the comment above)
!            [a],[b],.. element matrices needed to compute elk
!            [elk]..... element coefficient matrix [k]
!            [elm]..... element ’mass’ matrix [m]
!         __________________________________________________________________

  !implicit real*8(a-h,o-z)

  include 'stf1.h'
  include 'stf2.h'
  include 'shp.h'

  real (kind=dbl) :: gauspt(5,5),gauswt(5,5),xi,x,sum,h,eij,dij,dji,di, &
                     denom,d22,d12,d11,cx,const,cij,c22,c12,c11,bx,bij, &
                     b11,b10,b01,b00,ax,anu21,a33,aij,fx=0.0d0,ct=0.0d0
  integer :: nn,ni,ngp,lgp,jj,j,i,ii

  data gauspt/5*0.0d0,-.57735027d0,.57735027d0,3*0.0d0,-.77459667d0, &
  0.0d0,.77459667d0,2*0.0d0,-.86113631d0,-.33998104d0,.33998104d0,   &
  .86113631d0,0.0d0,-.906180d0,-.538469d0,0.0d0,.538469d0,.906180d0/

  data gauswt/2.0d0,4*0.0d0,2*1.0d0,3*0.0d0,.55555555d0,.88888888d0, &
  0.55555555d0,2*0.0d0,.34785485d0,2*.65214515d0,.34785485d0,0.0d0,  &
  0.236927d0,.478629d0,.568889d0,.478629d0,.236927d0/

  nn=ndf*npe
  h = elx(npe) - elx(1)
  if(ielem == 0)then
    ngp=4
  else
    ngp = ielem+1
  endif

  do j=1,nn
    if(item  <=  2)then
      elf(j) = 0.0
    endif
    do i=1,nn
      if(item >  0)then
        elm(i,j)=0.0
      endif
      elk(i,j)=0.0
    enddo
  enddo

  if(model /= 2)then

    ! do-loop on number of gauss points begins here
    do ni=1,ngp  ! do 100
      xi = gauspt(ni,ngp)
      ! call subroutine shape1d to evaluate the interpolation functions
      ! and their global derivatives at the gauss point xi

      call shape1d(h,ielem,npe,xi)
      const = gj*gauswt(ni,ngp)

      if(ielem == 0)then
        x = elx(1) + 0.5*h*(1.0+xi)
      else
        x = 0.0
        do j=1,npe
          x = x + sf(j)*elx(j)
        enddo
      endif

      ! compute coefficient matrices and vectors for vaious model problems
      ! governed by single second-order and fourth-order equations
      ! (model = 1 or 3; ntype = 0 or 1)

      cx=cx0+cx1*x
      if(item /= 3) then
         fx=fx0+fx1*x+fx2*x*x
      endif
      if(item > 0)then
         ct=ct0+ct1*x
      endif
      if(model == 1)then
        ! coefficients for all single-variable problems (model=1)
        if(ntype == 0)then
          ! all problems governed by model equation (3.1) (ntype=0)
          ax=ax0+ax1*x
          do j = 1,nn
            if(item <= 2)then
              elf(j) = elf(j) + const*sf(j)*fx
            endif
            do i = 1,nn
              if(item /= 0)then
                elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct
              endif
              aij = const*gdsf(i)*gdsf(j)
              cij = const*sf(i)*sf(j)
              elk(i,j)=elk(i,j) + ax*aij + cx*cij
            enddo
          enddo
        else
          ! radially symmetric elasticity problems (model=1, ntype>0)
          ! ax0=e1, ax1=e2, bx0=nu12, bx1=h, thickness

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


          do j=1,nn
            if(item <= 2)then
              elf(j) = elf(j) + const*sf(j)*fx*x
            endif
            do i=1,nn
              if(item /= 0)then
                elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct*x
              endif
              aij = const*gdsf(i)*gdsf(j)*c11*x
              cij = const*sf(i)*sf(j)*cx*x
              dij = const*(gdsf(i)*sf(j)+sf(i)*gdsf(j))*c12
              eij = const*sf(i)*sf(j)*c22/x
              elk(i,j)=elk(i,j) + aij + cij + dij + eij
            enddo
          enddo
        endif
      else
        ! coefficients for the euler-bernoulli theory (model=2)

        if(ntype == 0)then
          ! the euler-bernoulli beam element (model=1 and ntype=0)
          bx=bx0+bx1*x
          cx=cx0+cx1*x
          do j = 1,nn
            if(item <= 2)then
              elf(j) = elf(j) + const*sf(j)*fx
            endif
            do i = 1,nn
              if(item > 0)then
                if(item <= 3)then
                  elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct
                else
                  elm(i,j) = elm(i,j) + const*gdsf(i)*gdsf(j)
                endif
              endif
              bij = const*gddsf(i)*gddsf(j)
              cij = const*sf(i)*sf(j)
              elk(i,j)=elk(i,j) + bx*bij + cx*cij
            enddo
          enddo
        else
          ! the e-b circular plate element (model=1 and ntype>0)
          anu21=bx0*ax0/ax1
          di=(bx1**3)/12.0
          d11=di*ax0/(1.0-bx0*anu21)
          d22=d11*(ax1/ax0)
          d12=bx0*d22
          do j=1,nn
            if(item <= 2)then
              elf(j) = elf(j) + const*sf(j)*fx*x
            endif
            do i=1,nn
              bij = const*gddsf(i)*gddsf(j)*d11*x
              cij = const*sf(i)*sf(j)*cx*x
              dij = const*(gddsf(i)*gdsf(j)+gdsf(i)*gddsf(j))*d12
              eij = const*gdsf(i)*gdsf(j)*d22/x
              elk(i,j)=elk(i,j) + bij + cij + dij + eij
            enddo
          enddo
        endif
      endif
    enddo ! 100
  else  ! this is correct
    ! coefficients for the timoshenko beam and circular plate (model=2)
    ! full integration for bending coefficients
    do ni=1,ngp
      xi=gauspt(ni,ngp)
      call shape1d(h,ielem,npe,xi)
      const=gj*gauswt(ni,ngp)
      x = 0.0
      do j=1,npe
        x = x + sf(j)*elx(j)
      enddo

      if(ntype == 0 .or. ntype == 2)then
        ! the timoshenko beam element (model=2 and ntype=0 or 2)
        bx=bx0+bx1*x
        cx=cx0+cx1*x
        fx=fx0+fx1*x+fx2*x*x
        jj=1
        do j=1,npe
          if(item <= 2)then
            elf(jj)=elf(jj)+fx*sf(j)*const
          endif
          ii=1
          do  i=1,npe
            cij=sf(i)*sf(j)*const
            bij=gdsf(i)*gdsf(j)*const
            elk(ii,jj)    =elk(ii,jj)    +cx*cij
            elk(ii+1,jj+1)=elk(ii+1,jj+1)+bx*bij
            if(item /= 0)then
              elm(ii,jj)    =elm(ii,jj)    +ct0*cij
              elm(ii+1,jj+1)=elm(ii+1,jj+1)+ct1*cij
            endif
            ii=ndf*i+1
          enddo
          jj=ndf*j+1
        enddo
      else

        ! timoshenko circular plate element (model=2 and ntype=1 or 3)
        ! ax0=e1, ax1=e2, bx0=anu12, bx1=h

        anu21=bx0*ax0/ax1
        cx=cx0+cx1*x
        fx=fx0+fx1*x
        di=(bx1**3)/12.0
        d11=di*ax0/(1.0-bx0*anu21)
        d22=d11*(ax1/ax0)
        d12=bx0*d22
        jj=1
        do j=1,npe
          if(item <= 2)then
            elf(jj)=elf(jj)+fx*sf(j)*const*x
          endif
          ii=1
          do i=1,npe
            bij = const*gdsf(i)*gdsf(j)*d11*x
            cij = const*sf(i)*sf(j)*x
            dij = const*(gdsf(i)*sf(j)+sf(i)*gdsf(j))*d12
            eij = const*sf(i)*sf(j)*d22/x
            elk(ii,jj)    =elk(ii,jj)     + cx*cij
            elk(ii+1,jj+1)=elk(ii+1,jj+1) + bij + dij + eij
            if(item /= 0)then
              elm(ii,jj)    =elm(ii,jj)    +ct0*cij
              elm(ii+1,jj+1)=elm(ii+1,jj+1)+ct1*cij
            endif
            ii=ndf*i+1
          enddo
          jj=ndf*j+1
        enddo
      endif
    enddo

    ! reduced integration is used to evaluate the transverse shear terms
    lgp=ngp-1
    do ni=1,lgp
      xi=gauspt(ni,lgp)

      call shape1d(h,ielem,npe,xi)
      const=gj*gauswt(ni,lgp)

      x = 0.0
      do j=1,npe
        x = x + sf(j)*elx(j)
      enddo
      if(ntype == 0 .or. ntype == 2)then
        ! the timoshenko beam element (model=2 and ntype=0 or 2)
        ! ax = gak = ax0 + ax1*x (reduced integration)

        ax=ax0+ax1*x
        jj=1
        do j=1,npe
          ii=1
          do i=1,npe
            b11=gdsf(i)*gdsf(j)*const
            b01=sf(i)*gdsf(j)*const
            b10=gdsf(i)*sf(j)*const
            b00=sf(i)*sf(j)*const
            elk(ii,jj)    =elk(ii,jj)    +ax*b11
            elk(ii,jj+1) =elk(ii,jj+1) +ax*b10
            elk(ii+1,jj) =elk(ii+1,jj) +ax*b01
            elk(ii+1,jj+1)=elk(ii+1,jj+1)+ax*b00
            ii=i*ndf+1
          enddo
          jj=j*ndf+1
        enddo
      else
        ! timoshenko circular plate element (model=2 and ntype=1 or 3)
        ! bx1=h, fx2=g13*k (reduced integration)
        a33=bx1*fx2
        jj=1
        do j=1,npe
          ii=1
          do i=1,npe
            bij = const*gdsf(i)*gdsf(j)*x
            cij = const*sf(i)*sf(j)*x
            dij = const*gdsf(i)*sf(j)*x
            dji = const*sf(i)*gdsf(j)*x
            elk(ii,jj)    =elk(ii,jj)     + a33*bij
            elk(ii,jj+1) =elk(ii,jj+1)    + a33*dij
            elk(ii+1,jj) =elk(ii+1,jj)    + a33*dji
            elk(ii+1,jj+1)=elk(ii+1,jj+1) + a33*cij
            ii=ndf*i+1
          enddo
          jj=ndf*j+1
        enddo
      endif
    enddo

    if(item == 0 .and. ntype > 1)then
      call timforce(elf,elx,fx0,fx1,fx2,h,ntype,ne,f3,mxelm)
    endif
  endif

  if(item > 2)return

  if(item == 1 .or. item == 2)then
    ! equivalent coefficient matrices for time-dependent problems
    if(item  ==  1)then
      ! alfa-family of time approximation for parabolic equations
      do j=1,nn
        sum=0.0
        do i=1,nn
          sum=sum+(elm(i,j)-a2*elk(i,j))*elu(i)
          elk(i,j)=elm(i,j)+a1*elk(i,j)
        enddo
        elf(j)=(a1+a2)*elf(j)+sum
      enddo
    else
      ! newmark-family of time approximation for hyperbolic equations
      if(time == 0.0)then
        do j=1,nn
          do i=1,nn
            elf(j)=elf(j)-elk(i,j)*elu(i)
            elk(i,j)=elm(i,j)
          enddo
        enddo
      else
        do j=1,nn
          do i=1,nn
            elf(j)=elf(j)+elm(i,j)*(a3*elu(i)+a4*elv(i)+a5*ela(i))
            elk(i,j)=elk(i,j)+a3*elm(i,j)
          enddo
        enddo
      endif
    endif
  endif
  return
end ! end of subroutine coeffcnt() !
