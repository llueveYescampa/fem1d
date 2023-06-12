subroutine boundary(neq,neqr,nhbw,nspv,nssv,nnbc,ndf,dt,item,alfa,ibdy,ispv,issv, &
                    inbc,uref,vspv,vssv,vnbc,glk,glm,glf,gu0,mxebc,mxnbc,mxmbc,mxneq)
include 'common.h'

integer, intent(in)           :: neq,nhbw,nspv,nssv,nnbc,ndf, item, mxebc,mxnbc,mxmbc,mxneq
integer, intent(in)           :: ispv(mxebc,2),issv(mxnbc,2),inbc(mxmbc,2)
integer, intent(inout)        :: neqr, ibdy(mxebc)

real (kind=dbl) ,intent(in)   :: dt, alfa, uref(mxmbc),vspv(mxebc),vssv(mxnbc),vnbc(mxmbc), gu0(mxneq)
real (kind=dbl) ,intent(inout):: glk(mxneq,mxneq), glm(mxneq,mxneq), glf(mxneq)


!__________________________________________________________________
!
!   the subroutine is called in main to implement specified boundary
!   conditions on the assembled system of finite element equations
!__________________________________________________________________

  integer         :: ib,nf,neqr1,nb,jj,j,i,ic,ie,ii,ikept,imax,it,nc

!           impose boundary conditions for static and time-dependent problems

  if(item  <= 2) then
!   include specified primary degrees of freedom

    if(nspv /=  0) then
      do nb = 1,nspv
        ie=(ispv(nb,1)-1)*ndf+ispv(nb,2)
        it=nhbw-1
        i=ie-nhbw
        do ii=1,it
          i=i+1
          if(i >= 1) then
            j=ie-i+1
            glf(i)=glf(i)-glk(i,j)*vspv(nb)
            glk(i,j)=0.0
          endif
        enddo
        glk(ie,1)=1.0
        glf(ie)=vspv(nb)
        i=ie
        do ii=2,nhbw
          i=i+1
          if(i  <=  neq)then
            glf(i)=glf(i)-glk(ie,ii)*vspv(nb)
            glk(ie,ii)=0.0
          endif
        enddo
      enddo
    endif

    if (nssv  /=  0) then
      ! include specified secondary degrees of freedom
      do nf = 1,nssv
        nb=(issv(nf,1)-1)*ndf+issv(nf,2)
        if(item  == 1) then
          glf(nb)=glf(nb)+vssv(nf)*dt
        else
          glf(nb)=glf(nb)+vssv(nf)
        endif
      enddo
    endif

    if (nnbc /=  0) then
    ! include specified mixed boundary conditions
      do ic=1,nnbc
        nc=(inbc(ic,1)-1)*ndf+inbc(ic,2)
        if(item == 1)then
          glk(nc,1)=glk(nc,1)+alfa*dt*vnbc(ic)
          glf(nc)=glf(nc)+dt*vnbc(ic)*(uref(ic)-(1.0-alfa)*gu0(nc))
        else
          glk(nc,1)=glk(nc,1)+vnbc(ic)
          glf(nc)=glf(nc)+vnbc(ic)*uref(ic)
        endif
      enddo
    endif

  else

! impose boundary conditions for eigenvalue problems

    if(nnbc /= 0) then
    !  include specified mixed boundary conditions
      do ic=1,nnbc
        nc=(inbc(ic,1)-1)*ndf+inbc(ic,2)
        glk(nc,nc)=glk(nc,nc)+vnbc(ic)
      enddo
    endif

!   include specified primary degrees of freedom

    if(nspv /=  0)then
      do  ib=1,nspv
        ibdy(ib)=(ispv(ib,1)-1)*ndf+ispv(ib,2)
      enddo
      do i=1,nspv
        imax=ibdy(i)
        do j=i,nspv
          if(ibdy(j).ge.imax)then
            imax=ibdy(j)
            ikept=j
          endif
        enddo
        ibdy(ikept)=ibdy(i)
        ibdy(i)=imax
      enddo
      neqr = neq
      do i=1,nspv
        ib=ibdy(i)
        if(ib  <  neqr)then
          neqr1=neqr-1
          do ii=ib,neqr1
            do jj=1,neqr
              glm(ii,jj)=glm(ii+1,jj)
              glk(ii,jj)=glk(ii+1,jj)
            enddo
            do jj=1,neqr
              glm(jj,ii)=glm(jj,ii+1)
              glk(jj,ii)=glk(jj,ii+1)
            enddo
          enddo
        endif
        neqr=neqr-1
      enddo
    endif
  endif
  return
end ! end of subroutine boundary() !
