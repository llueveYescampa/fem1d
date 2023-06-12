subroutine transfrm(mxelm,n,ntype,pr,se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
include 'common.h'

integer, intent(in)         :: mxelm,n,ntype
real (kind=dbl) ,intent(in) :: pr(mxelm),se(mxelm),sl(mxelm),sa(mxelm),si(mxelm),cs(mxelm),sn(mxelm), &
                               cnt(mxelm),snt(mxelm),hf(mxelm),vf(mxelm),pf(mxelm),xb(mxelm)


!        __________________________________________________________________

!           called in both main and reaction to compute stiffness matrix and
!            force vector for the truss (ndf=2) and frame (ndf=3) elements

!         se......young’s modulus
!         sl......element length
!         sa......cross-sectional area
!         si......moment of inertia
!         cs......cosine of the angle of orientation
!         sn......sine of the angle of orientation
!         hf......distributed force along the length of the element
!         vf......distributed force transverse to the element
!         pf......point force at point other than nodes
!         xb......distance along the length from node 1 of the element
!                 of the location of the point force, pf
!         cnt,snt:direction cosines of the point force’s line of application
!         __________________________________________________________________


  real (kind=dbl) :: trm(6,6),tmpk(6,6),c1,c2,c3,c4,c5,c6,cn1,cn2,csn,f1,f2,f3,f4,f5,f6, &
                     af,amu,sfh1,sfh2, sfh3,sfh4,sfl1,sfl2,sfq1,sfq2,sfq3,sg,sn1,sn2,tf,xi
  integer         :: i,j,k,nn

  include 'stf1.h'
  include 'io.h'

  cn1=cs(n)
  sn1=sn(n)
  cn2=cn1*cn1
  sn2=sn1*sn1
  csn=cn1*sn1
  nn = 6 ! dof for case of frame elements used as default
!element coefficients


  select case (ntype)
  case (0)
    ! the plane truss element

    nn=4 ! dof for case of plane truss element
    c1=sa(n)*se(n)/sl(n)
    elk(1,1) = c1*cn2
    elk(2,1) = c1*csn
    elk(2,2) = c1*sn2
    elk(3,1) = -elk(1,1)
    elk(3,2) = -elk(2,1)
    elk(3,3) = elk(1,1)
    elk(4,1) = -elk(2,1)
    elk(4,2) = -elk(2,2)
    elk(4,3) = -elk(3,2)
    elk(4,4) = elk(2,2)

    do i=1,nn
      do j=i,nn
        elk(i,j) = elk(j,i)
      enddo
    enddo

  ! contribution of the point force to nodal forces

    xi=xb(n)/sl(n)
    sfl1 = 1.0-xi
    sfl2 = xi

    f1=0.5*hf(n)*sl(n)
    f3=0.5*hf(n)*sl(n)
    elf(1) = f1*cn1
    elf(2) = f1*sn1
    elf(3) = f3*cn1
    elf(4) = f3*sn1

  case (1)
    ! the euler-bernoulli frame element
    amu=0.5*sa(n)*sl(n)*sl(n)/si(n)
    c1=2.0*se(n)*si(n)/(sl(n)**3)
    c2=6.0*se(n)*si(n)/(sl(n)*sl(n))
    c3=c1*(amu*cn2+6.0*sn2)
    c4=c1*(amu-6.0)*csn
    c5=c1*(amu*sn2+6.0*cn2)
    c6=4.0*se(n)*si(n)/sl(n)

    elk(1,1)   = c3
    elk(2,1)   = c4
    elk(2,2)   = c5
    elk(3,1)   = c2*sn1
    elk(3,2)   =-c2*cn1
    elk(3,3)   = c6
    elk(4,1)   =-c3
    elk(4,2)   =-c4
    elk(4,3)   =-c2*sn1
    elk(4,4)   = c3
    elk(5,1)   =-c4
    elk(5,2)   =-c5
    elk(5,3)   = c2*cn1
    elk(5,4)   = c4
    elk(5,5)   = c5
    elk(6,1)   = c2*sn1
    elk(6,2)   =-c2*cn1
    elk(6,3)   = 0.5*c6
    elk(6,4)   =-c2*sn1
    elk(6,5)   = c2*cn1
    elk(6,6)   = c6

    do i=1,nn
      do j=i,nn
        elk(i,j) = elk(j,i)
      enddo
    enddo

    ! contribution of the point force to nodal generalized forces

    xi=xb(n)/sl(n)
    tf=pf(n)*snt(n)
    af=pf(n)*cnt(n)
    sfl1 = 1.0-xi
    sfl2 = xi
    sfh1 = 1.0 - 3.0*xi*xi + 2.0*(xi**3)
    sfh2 = -xi*(1.0+xi*xi-2.0*xi)*sl(n)
    sfh3 = 3.0*xi*xi - 2.0*(xi**3)
    sfh4 = -xi*(xi*xi - xi)*sl(n)

    f1=0.5*hf(n)*sl(n)                     +   sfl1*af
    f2=0.5*vf(n)*sl(n)                     +   sfh1*tf
    f3=-vf(n)*sl(n)*sl(n)/12.0             +   sfh2*tf
    f4=0.5*hf(n)*sl(n)                     +   sfl2*af
    f5=0.5*vf(n)*sl(n)                     +   sfh3*tf
    f6=vf(n)*sl(n)*sl(n)/12.0              +   sfh4*tf
    elf(1) = f1*cn1-f2*sn1
    elf(2) = f1*sn1+f2*cn1
    elf(3) = f3
    elf(4) = f4*cn1-f5*sn1
    elf(5) = f4*sn1+f5*cn1
    elf(6) = f6

  case (2)
    ! the timoshenko frame element (shear coefficient=5/6)

    sg=5.0*se(n)/(1.0+pr(n))/12.0
    c1=sa(n)*se(n)/sl(n)
    c2=sg*sa(n)/sl(n)
    c3=0.5*sg*sa(n)
    c4=0.25*sg*sa(n)*sl(n)
    c5=se(n)*si(n)/sl(n)
    elk(1,1)=c1
    elk(2,1)=0.0
    elk(2,2)=c2
    elk(3,1)=0.0
    elk(3,2)=-c3
    elk(3,3)=c4+c5
    elk(4,1)=-c1
    elk(4,2)=0.0
    elk(4,3)=0.0
    elk(4,4)=c1
    elk(5,1)=0.0
    elk(5,2)=-c2
    elk(5,3)=c3
    elk(5,4)=0.0
    elk(5,5)=c2
    elk(6,1)=0.0
    elk(6,2)=-c3
    elk(6,3)=c4-c5
    elk(6,4)=0.0
    elk(6,5)=c3
    elk(6,6)=c4+c5

    do i=1,nn
      do j=1,nn
        trm(j,i)=0.0
      enddo
    enddo

    trm(1,1)=cn1
    trm(1,2)=sn1
    trm(2,1)=-sn1
    trm(2,2)=cn1
    trm(3,3)=1.0
    trm(4,4)=cn1
    trm(4,5)=sn1
    trm(5,4)=-sn1
    trm(5,5)=cn1
    trm(6,6)=1.0


    do i=1,nn
      do j=i,nn
        elk(i,j) = elk(j,i)
      enddo
    enddo

    do i=1,nn
      do j=1,nn
        tmpk(i,j)=0.0
        do k=1,nn
          tmpk(i,j)=tmpk(i,j)+trm(k,i)*elk(k,j)
        enddo
      enddo
    enddo

    do i=1,nn
      do j=1,nn
        elk(i,j)=0.0
        do k=1,nn
          elk(i,j)=elk(i,j)+tmpk(i,k)*trm(k,j)
        enddo
        enddo
    enddo

    ! contribution of the point force to nodal generalized forces

    xi=xb(n)/sl(n)
    tf=pf(n)*snt(n)
    af=pf(n)*cnt(n)
    sfl1 = 1.0-xi
    sfl2 = xi
    sfq1 = (1.0-xi)*(1.0-2.0*xi)
    sfq2 = -xi*(1.0-2.0*xi)
    sfq3 = 4.0*xi*(1.0-xi)

    f1=0.5*hf(n)*sl(n)                     +   sfl1*af
    f2=0.5*vf(n)*sl(n)                     +   (sfq1+0.5*sfq3)*tf
    f3=-vf(n)*sl(n)*sl(n)/12.0             -   0.125*sfq3*sl(n)*tf
    f4=0.5*hf(n)*sl(n)                     +   sfl2*af
    f5=0.5*vf(n)*sl(n)                     +   (sfq2+0.5*sfq3)*tf
    f6=vf(n)*sl(n)*sl(n)/12.0              +   0.125*sfq3*sl(n)*tf
    elf(1) = f1*cn1-f2*sn1
    elf(2) = f1*sn1+f2*cn1
    elf(3) = f3
    elf(4) = f4*cn1-f5*sn1
    elf(5) = f4*sn1+f5*cn1
    elf(6) = f6
  case default
    write(it, "('Element type not recognized. Stoping ....',/)"  )
    stop
  end select
  return
end ! end of subroutine transfrm() !
