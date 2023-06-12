subroutine jacobi (n,q,jvec,m,v,x,ih,mxneq)
include 'common.h'

integer, intent(in)            :: n,jvec,mxneq
integer, intent(inout)         :: m, ih(mxneq)
real (kind=dbl) ,intent(inout) :: x(mxneq),v(mxneq,mxneq), q(mxneq,mxneq)

!        __________________________________________________________________

!          called in egnsolvr to diagonalize [q] by successive rotations

!          description of the variables:

!            n   .... order of the real, symmetric matrix [q] (n > 2)
!           [q] .... the matrix to be diagonalized (destroyed)
!           jvec .... 0, when only eigenvalues alone have to be found
!           [v] .... matrix of eigenvectors
!            m   .... number of rotations performed
!        __________________________________________________________________

  !implicit real*8 (a-h,o-z)

  integer         :: mj,mi,i,j,k
  integer         :: ip=0, jp=0
  real (kind=dbl) :: temp, sine, tang,qii,xmax=0.0,cosn
  real (kind=dbl), parameter :: epsi=1.0d-08


  if (jvec /= 0) then
    do i=1,n
      do j=1,n
        if (i == j) then
          v(i,j)=1.0
        else
          v(i,j)=0.0
        endif
      enddo
    enddo
  endif

  m=0
  mi=n-1

  do i=1,mi
    x(i)=0.0
    mj=i+1
    do j=mj,n
      if ( x(i) <= dabs(q(i,j)) ) then
        x(i)=dabs(q(i,j))
        ih(i)=j
      endif
    enddo
  enddo

  do ! while
    do i=1,mi
      if (i <= 1 .or. (i> 1 .and. xmax < x(i)) ) then
        xmax=x(i)
        ip=i
        jp=ih(i)
      end if
    enddo
    if (xmax>epsi) then
      m=m+1
      if ( q(ip,ip) < q(jp,jp) ) then
        tang=-2.0*q(ip,jp)/(dabs(q(ip,ip)-q(jp,jp))+dsqrt((q(ip,ip) - q(jp,jp))**2 + 4.0*q(ip,jp)**2))
      else
        tang= 2.0*q(ip,jp)/(dabs(q(ip,ip)-q(jp,jp))+dsqrt((q(ip,ip) - q(jp,jp))**2 + 4.0*q(ip,jp)**2))
      endif
      cosn=1.0/dsqrt(1.0+tang**2)
      sine=tang*cosn
      qii=q(ip,ip)
      q(ip,ip)=cosn**2*(qii+tang*(2.*q(ip,jp)+tang*q(jp,jp)))
      q(jp,jp)=cosn**2*(q(jp,jp)-tang*(2.*q(ip,jp)-tang*qii))
      q(ip,jp)=0.0
      if ( q(ip,ip) <  q(jp,jp)  ) then
        temp=q(ip,ip)
        q(ip,ip)=q(jp,jp)
        q(jp,jp)=temp
        if ( sine < 0.0 ) then
          temp=cosn
        else
          temp=-cosn
        endif
        cosn=dabs(sine)
        sine=temp
      endif
      do i=1,mi
         if ( ( i > ip .and. i/=jp .and. ( ( ih(i) == ip .or. ih(i) /= ip .and. ih(i) == jp )  ) )  .or. &
              ( i < ip .and.             ( ( ih(i) == ip .or. ih(i) /= ip .and. ih(i) == jp )  ) ) ) then
           k=ih(i)
           temp=q(i,k)
           q(i,k)=0.0
           mj=i+1
           x(i)=0.0
           do j=mj,n
             if ( x(i)  <=  dabs(q(i,j))  ) then
               x(i)=dabs(q(i,j))
               ih(i)=j
             endif
           enddo
           q(i,k)=temp
         endif
      enddo

      x(ip)=0.0
      x(jp)=0.0
      do i=1,n
        if (i < ip) then
          temp=q(i,ip)
          q(i,ip)=cosn*temp+sine*q(i,jp)
          if ( x(i) < dabs(q(i,ip)) ) then
            x(i)=dabs(q(i,ip))
            ih(i)=ip
          endif
          q(i,jp)=-sine*temp+cosn*q(i,jp)
          if ( x(i) < dabs(q(i,jp)) ) then
            x(i)=dabs(q(i,jp))
            ih(i)=jp
          endif
        else if (i > ip) then
          if (i < jp  ) then
            temp=q(ip,i)
            q(ip,i)=cosn*temp+sine*q(i,jp)
            if (x(ip) < dabs(q(ip,i)) ) then
              x(ip)=dabs(q(ip,i))
              ih(ip)=i
            endif
            q(i,jp)=-sine*temp+cosn*q(i,jp)

            if (x(i) < dabs(q(i,jp))) then
              x(i)=dabs(q(i,jp))
              ih(i)=jp
            endif
          else if ( i > jp ) then
            temp=q(ip,i)
            q(ip,i)=cosn*temp+sine*q(jp,i)
            if(x(ip) < dabs(q(ip,i))) then
              x(ip)=dabs(q(ip,i))
              ih(ip)=i
            endif
            q(jp,i)=-sine*temp+cosn*q(jp,i)
            if(x(jp) < dabs(q(jp,i))) then
              x(jp)=dabs(q(jp,i))
              ih(jp)=i
            end if
          endif
        endif
      enddo

      if (jvec /= 0) then
        do i=1,n
          temp=v(i,ip)
          v(i,ip)=cosn*temp+sine*v(i,jp)
          v(i,jp)=-sine*temp+cosn*v(i,jp)
        enddo
      endif
    else
      exit
    endif
  enddo
  return
end ! end of subroutine jacobi() !
