subroutine reaction(mxelm,mxneq,ndf,nem,nod,npe,ntype,pr,glf,se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
include 'common.h'

integer, intent(in)           :: mxelm,mxneq,ndf,nem,nod(mxelm,4),npe,ntype
real (kind=dbl) ,intent(in)   :: pr(mxelm),glf(mxneq),se(mxelm),sl(mxelm),sa(mxelm),si(mxelm),xb(mxelm)
real (kind=dbl) ,intent(in)   :: cs(mxelm),sn(mxelm),cnt(mxelm),snt(mxelm),hf(mxelm),vf(mxelm),pf(mxelm)

!        __________________________________________________________________

!         the subroutine is called in main to compute generalized reaction
!        forces in each element of truss (ndf=2) or frame (ndf=3) structure
!        __________________________________________________________________

  real (kind=dbl) ::elr(6),cn1,sn1
  integer         :: n,l,nn,li,i,j,ni

  include 'stf1.h'

  nn=npe*ndf
  do n=1,nem
    cn1=cs(n)
    sn1=sn(n)

    ! call transfrm to compute element stiffness matrix and force vector

    l=0
    do i=1,npe
      ni=nod(n,i)
      li=(ni-1)*ndf
      do j=1,ndf
        li=li+1
        l=l+1
        elu(l)=glf(li)
      enddo
    enddo
    call transfrm(mxelm,n,ntype,pr,se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)

    !  compute the force and moment resultants

    do i=1,nn
      elr(i) = 0.0
      do j=1,nn
        elr(i) = elr(i) + elk(i,j)*elu(j)
      enddo
      elr(i) = elr(i) - elf(i)
    enddo

    elf(1) = elr(1)*cn1+elr(2)*sn1
    elf(2) = -elr(1)*sn1+elr(2)*cn1
    if(ntype.ne.0) then
      elf(3) = elr(3)
      elf(4) = elr(4)*cn1+elr(5)*sn1
      elf(5) = -elr(4)*sn1+elr(5)*cn1
      elf(6) = elr(6)
    else
      elf(3) = elr(3)*cn1+elr(4)*sn1
      elf(4) = -elr(3)*sn1+elr(4)*cn1
    endif
    write(6, "(3x,i2,6e12.4)" ) n, (elf(i),i=1,nn)
    write(6, "(5x,6e12.4,/)")   (elr(i),i=1,nn)

  enddo
  return
end ! end of subroutine reaction() !
