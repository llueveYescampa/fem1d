!       Program Name: FEM1D           Length(INCLUDINMG BLANKS):2440 lines
!
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       *                        Program FEM1D                          *
!       *          (A FINITE ELEMENT ANALYSIS COMPUTER PROGRAM)         *
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         _______________________________________________________________
!       |                                                                |
!       | This is a finite element computer program for the analysis     |
!       | of the following three model equations and others:             |
!       |                                                                |
!       | 1. Heat transfer, fluid mechanics, bars, and cables:           |
!       |                                                                |
!       |             CT.u* + CT.u** - (AX.u')' + CX.u = FX              |
!       |                                                                |
!       | 2. The Timoshenko beam and circular plate theory:              |
!       |                                                                |
!       |              CT0.w** - [AX.(w' + s)]' + CX.w = FX              |
!       |              CT1.s** - (BX.s')' + AX.(w' + s) = 0              |
!       |                                                                |
!       | 3. The Euler-Bernoulli beam and circular plate theory:         |
!       |                                                                |
!       |                CT.w** + (BX.w'')'' + CX.w = FX                 |
!       |                                                                |
!       | In the above equations (') and (*) denote differentiations     |
!       | with respect to space x and time t, and AX, BX, CX, CT, and    |
!       | FX are functions of x only:                                    |
!       |                                                                |
!       |     AX = AX0 + AX1.X, BX = BX0 + BX1.X, CX = CX0 + CX1.X       |
!       |          CT = CT0 + CT1.X, FX = FX0 + FX1.X + FX2.X.X          |
!       |                                                                |
!       |     In addition to the three model equations, other equations  |
!       | (for example, disks, trusses, and frames) can be analyzed by   |
!       | the program.                                                   |
!       |________________________________________________________________|
!
!        _________________________________________________________________
!        .                                                               .
!        .               KEY VARIABLES USED IN THE PROGRAM               .
!        . See Table 7.3.2 of the BOOK for a description of the variables.
!        .                                                               .
!        . NDF....... Number of degrees of freedom per node              .
!        . NEQ....... Number of equations in the model (before B. C.)    .
!        . NGP....... Number of Gauss points used in the evaluation of   .
!        .            the element coefficients, ELK , ELF , ELM          .
!        . NHBW...... Half bandwidth of global coefficient matrix GLK    .
!        . NN ....... Number of total degrees of freedom in the element .
!        . NPE....... Number of nodes per element                        .
!        _________________________________________________________________
!        _________________________________________________________________
!        .                                                               .
!        .          DIMENSIONS OF VARIOUS ARRAYS IN THE PROGRAM          .
!        .                                                               .
!        . Values of MXELM,MXNOD, etc. in the PARAMETER statement should .
!        .         be changed to meet the requirements of the problem:   .
!        .                                                               .
!        . MXELM..... Maximum number of elements in the mesh:            .
!        . MXEBC..... Maximum number of speci. primary deg. of freedom   .
!        . MXMBC..... Maximum number of speci. mixed boundary conditions .
!        . MXNBC..... Maximum number of speci. secondary deg. of freedom .
!        . MXNEQ..... Maximum number of equations in the FE model        .
!        . MXNOD..... Maximum number of nodes in the mesh                .
!        .                                                               .
!        . NOTE: The following dimension statement in subroutine JACOBI .
!        .        should be modified when MXNEQ is greater than 500:     .
!        .         DIMENSION V(500,500),VT(500,500),W(500,500),IH(500) .
!        .        The value of MXNEQ should be used in place of ‘500'    .
!        _________________________________________________________________
!        .                                                               .
!        .                SUBROUTINES USED IN THE PROGRAM                .
!        .                                                               .
!        . ASSEMBLE, BOUNDARY, COEFFCNT, CONSTRNT, ECHODATA, EQNSOLVR, .
!        .    EIGNSLVR, JACOBI, MATRXMLT, MESH1D, POSTPROC, REACTION,    .
!        .             SHAPE1D, TIMFORCE, TIMSTRES, TRANSFRM             .
!        _________________________________________________________________
!
program fem1d
Use, intrinsic :: iso_fortran_env, Only : iostat_end
include 'common.h'

!implicit real*8(a-h,o-z)
!implicit none
  integer           :: argc
  character(len=32) :: argv

  include 'formats4fem1d.f95'

  integer, parameter :: mxelm=250,mxneq=500,mxebc=20,mxnbc=20,mxmbc=20,mxnod=250,mxmpc=5

  !!!!!!!!!!!!!!!!!!!! Dynamically allocated variables !!!!!!!!!!!!!!!!!!!
  integer, dimension(:,:), allocatable :: ispv,issv,inbc,imc1,imc2,nod
  integer, dimension(:),   allocatable :: ibdy,icon

  real (kind=dbl), dimension(:,:), allocatable :: dcax,dcbx,dccx,dcfx,glm,trm,vmpc,egnvec,glk
  real (kind=dbl), dimension(:),   allocatable :: gu0,gu1,gu2,gpu,dx,vcon, glf,glx,pr,se,sl,vssv,vnbc
  real (kind=dbl), dimension(:),   allocatable :: cs,sn,cnt,snt,xb,hf,vf,pf,f3,sa,si,egnval,uref,vspv
  !!!!!!!!!!!!!!!!!!!! Dynamically allocated variables !!!!!!!!!!!!!!!!!!!

  integer            :: i,ib,icont,ielem,incond,intvl,ioerr,ires,item,j,jvec,l,li,model
  integer            :: n,nb,nc,ncon,ndf,ndof1,ndof2,ne,nem,nem1,neq,neqr,nhbw,ni,nmpc,nn
  integer            :: nnbc,nnm,npe,nprnt,nrot,nspv,nssv,nt,nten,ntime,ntype,nvec,nw
  real (kind=dbl) :: acclrn,alfa,diff,dt,frqncy,gama,prcnt,soln,time,tolrns,value1,vmax,pnlty=0.0
  character(len=80) :: title

  include 'io.h'
  include 'stf1.h'
  include 'stf2.h'

  !!!!!!!!!!!!!!!!!!!! Dynamically allocating variables !!!!!!!!!!!!!!!!!!!
  allocate ( ibdy(mxebc) )
  allocate ( icon(9) )

  allocate ( ispv(mxebc,2) )
  allocate ( issv(mxnbc,2) )
  allocate ( inbc(mxmbc,2) )

  allocate ( imc1(mxmpc,2) )
  allocate ( imc2(mxmpc,2) )
  allocate ( nod(mxelm,4) )

  allocate ( dcax(mxelm,2) )
  allocate ( dcbx(mxelm,2) )
  allocate ( dccx(mxelm,2) )
  allocate ( dcfx(mxelm,3) )

  allocate ( egnvec(mxneq,mxneq) )
  allocate ( glk(mxneq,mxneq) )

  allocate ( glm(mxneq,mxneq) )
  allocate ( trm(mxneq,mxneq) )
  allocate ( vmpc(mxmpc,4) )


  allocate ( gu0(mxneq) )
  allocate ( gu1(mxneq) )
  allocate ( gu2(mxneq) )
  allocate ( gpu(mxneq) )
  allocate ( dx(mxnod)  )
  allocate ( vcon(9) )
  allocate ( glf(mxneq) )
  allocate ( glx(mxnod) )
  allocate ( cs(mxelm) )
  allocate ( sn(mxelm) )
  allocate ( cnt(mxelm) )
  allocate ( snt(mxelm) )
  allocate ( xb(mxelm) )


  allocate ( hf(mxelm) )
  allocate ( vf(mxelm) )
  allocate ( pf(mxelm) )
  allocate ( f3(mxelm) )

  allocate ( pr(mxelm) )
  allocate ( se(mxelm) )
  allocate ( sl(mxelm) )
  allocate ( sa(mxelm) )
  allocate ( si(mxelm) )
  allocate ( egnval(mxneq) )

  allocate ( uref(mxmbc) )
  allocate ( vspv(mxebc) )
  allocate ( vssv(mxnbc) )
  allocate ( vnbc(mxmbc) )

  !!!!!!!!!!!!!!!!!!!! Dynamically allocating variables !!!!!!!!!!!!!!!!!!!

!                  _______________________________________________
!                |                                                |
!                |        p r e p r o c e s s o r   u n i t       |
!                |________________________________________________|

  in=5
  it=6

  argc = command_argument_count()
  if (argc == 0) then
    call get_command_argument(0, argv)
    print*, 'Use: ',  argv, 'filename  [outfile] '
    stop
  else if (argc >= 1) then
    call get_command_argument(1, argv)
    open(in, file = argv, status = 'old', iostat=ioerr)
    if (ioerr > 0) then
      print*, 'the ' ,  argv, 'input file was not found'
      print*, 'Bye ...'
      stop
    endif
    if (argc > 1) then
      call get_command_argument(2, argv)
      open(it, file = argv, status = 'unknown', iostat=ioerr)
    endif
  endif



  nt=0
  nssv=0
  jvec=1
  time=0.0d0
  tolrns=1.0d-06
  call echodata(in,it)

  read(in,'(80a)') title
  read(in,*) model,ntype,item
  read(in,*) ielem,nem
  read(in,*) icont,nprnt

  if(model >= 3)then
    npe=2
    if(model == 4 .and. ntype >= 1)then
      ndf=3
    else
      ndf=2
    endif
    if(model == 4 .and. ntype == 2)then
      ielem=1
    else
      ielem=0
    endif
  else
    if(model == 2)then
      ndf=2
      if(ntype > 1) ielem=1
    else
      ndf=1
    endif
    npe=ielem+1
  endif


  if(model /= 4)then
!   data input for bar-like and beam problems (model=1,2, and 3)
    if(icont /= 0)then
      nnm = nem*(npe-1)+1
      nem1=nem + 1
      read(in,*) (dx(i), i=1,nem1)
      call mesh1d(nem,npe,nod,mxelm,mxnod,dx,glx)
      read(in,*) ax0,ax1
      read(in,*) bx0,bx1
      read(in,*) cx0,cx1
      if(item < 3)then
        read(in,*) fx0,fx1,fx2
      endif
    else
    ! read glx, nod, and element-wise continuous coefficients [dc.x]
      read(in,*)nnm
      do n=1,nem
        read(in,*) (nod(n,i),i=1,npe), glx(n)
        read(in,*) (dcax(n,i),i=1,2)
        read(in,*) (dcbx(n,i),i=1,2)
        read(in,*) (dccx(n,i),i=1,2)
        read(in,*) (dcfx(n,i),i=1,3)
      enddo
    endif
  else
!   input data for plane truss or frame structures (model=4)
    read(in,*)nnm
    if(ntype /= 0)then
      do n=1,nem
        read(in,*) pr(n),se(n),sl(n),sa(n),si(n),cs(n),sn(n)
        read(in,*) hf(n),vf(n),pf(n),xb(n),cnt(n),snt(n)
        read(in,*) (nod(n,i),i=1,2)
      enddo
    else
      do n=1,nem
        read(in,*) se(n),sl(n),sa(n),cs(n),sn(n),hf(n)
        read(in,*) (nod(n,i),i=1,2)
      enddo
    endif
    read(in,*) ncon
    if(ncon /= 0)then
      do i=1, ncon
        read(in,*) icon(i),vcon(i)
      enddo
    endif
  endif

  neq=nnm*ndf

! read data on boundary conditions of three kinds: dirichlet (pv)
! neumann (sv), and newton's (mixed) types

  read(in,*) nspv
  if(nspv /= 0)then
    do nb=1,nspv
      if(item > 2)then
        read(in,*) (ispv(nb,j),j=1,2)
      else
        read(in,*) (ispv(nb,j),j=1,2),vspv(nb)
      endif
    enddo
  endif

  if(item <= 2)then
    read(in,*) nssv
    if(nssv /= 0)then
      do ib=1,nssv
        read(in,*) (issv(ib,j),j=1,2),vssv(ib)
      enddo
    endif
  endif

  read(in,*) nnbc
  if(nnbc /= 0)then
    do i=1, nnbc
      read(in,*) (inbc(i,j),j=1,2),vnbc(i),uref(i)
    enddo
  endif

! read data on multi-point constraints

  read(in,*) nmpc
  if(nmpc /= 0)then
    do i=1, nmpc
      read(in,*)(imc1(i,j),j=1,2),(imc2(i,j),j=1,2),(vmpc(i,j),j=1,4)
    enddo
  endif

  if(item  /= 0)then
!   input data here for time-dependent problems
    if(item <= 3)then
      read(in,*) ct0,ct1
    endif
    if(item <= 2)then
      read(in,*) dt,alfa,gama
      read(in,*) incond,ntime,intvl
      a1=alfa*dt
      a2=(1.0-alfa)*dt
      if(incond /= 0)then
        read(in,*) (gu0(i),i=1,neq)
      else
        do i=1,neq
          gu0(i)=0.0
        enddo
      endif
      if(item == 2)then
        a3=2.0/gama/(dt*dt)
        a4=a3*dt
        a5=1.0/gama-1.0
        if(incond /= 0)then
          read(in,*) (gu1(i),i=1,neq)
        else
          do i=1,neq
            gu1(i)=0.0
            gu2(i)=0.0
          enddo
        endif
      endif
    endif
  endif

! ----------------------------------------------------------------
!  e n d       o f      t h e       i n p u t          d a t a
! ----------------------------------------------------------------

! compute the half bandwidth of the coefficient matrix glk

  nhbw=0.0
  do n=1,nem
    do i=1,npe
      do j=1,npe
        nw=(iabs(nod(n,i)-nod(n,j))+1)*ndf
        if(nhbw < nw) nhbw=nw
      enddo
    enddo
  enddo



! ----------------------------------------------------------------
! p r i n t    t h e    i n p u t    d a t a
! ----------------------------------------------------------------

  write(it,"(2x,55('_'),/)")
  write(it,"(8x,'output from program    fem1d   by j n reddy')")
  write(it,"(2x,55('_'),/)")
  write(it,'(80a)') title
  write(it,320) model,ntype
  write(it,350) ielem,ndf,nem,neq,nhbw,nspv,nssv,nnbc,nmpc

  if(item /= 0)then
    if(item <= 2)then
      write(it,"(/,4x,'time-dependent (transient) analysis ',/)")
      write(it,390) ct0,ct1,alfa,gama,dt,ntime,intvl
      if(incond /= 0)then
        write(it,"(/,3x, 'initial conditions on the primary variables:',/)")
        write(it,"(2x,5e13.5)") (gu0(i),i=1,neq)
        if(item == 2)then
          write(it,"(/,3x, 'initial cond. on time der. of primary variables:',/)")
          write(it,"(2x,5e13.5)") (gu1(i),i=1,neq)
        endif
      endif
    else
      write(it,"(/,4x,'e i g e n v a l u e a n a l y s i s',/)")
      if(item <= 3)then
        write(it,400) ct0,ct1
      endif
    endif
  endif

  if(nspv /= 0)then
    write(it,"(/,3x, 'boundary information on primary variables:',/)")
    do ib=1,nspv
      if(item <= 2)then
        write(it,"(5x,2i5,2e15.5)") (ispv(ib,j),j=1,2),vspv(ib)
      else
        write(it,"(5x,2i5,2e15.5)") (ispv(ib,j),j=1,2)
      endif
    enddo
  endif

  if(nssv /= 0)then
    write(it,"(/,3x, 'boundary information on secondary variables:',/)")
    do ib=1,nssv
      write(it,"(5x,2i5,2e15.5)") (issv(ib,j),j=1,2),vssv(ib)
    enddo
  endif


  if(nnbc /= 0)then
    write(it,"(/,3x, 'boundary information on mixed boundary cond.:',/)")
    do i=1,nnbc
      write(it,"(5x,2i5,2e15.5)") (inbc(i,j),j=1,2),vnbc(i),uref(i)
    enddo
  endif

  if(nmpc /= 0)then
    write(it,"(/,3x, 'multi-point constraint information:',/)")
    do i=1, nmpc
      write(it,"(5x,2i5,2x,2i5,/,5x,4e15.5)")(imc1(i,j),j=1,2), &
                            (imc2(i,j),j=1,2),(vmpc(i,j),j=1,4)
    enddo
  endif

  if(model /= 4)then
    if(icont == 1)then
      write(it,"(/,3x,'global coordinates of the nodes, {glx}:',/)")
      write(it,"(2x,5e13.5)") (glx(i),i=1,nnm)
      write(it,"(/,3x,'coefficients of the differential equation:',/)")
      if(model /= 3)then
        write(it,440) ax0,ax1,bx0,bx1,cx0,cx1,fx0,fx1,fx2
      else
        write(it,445) ax0,ax1,bx0,bx1,cx0,cx1
      endif
    else
      do n=1,nem
        write(it,430) n,glx(n)
        write(it,440) (dcax(n,i),i=1,2),(dcbx(n,i),i=1,2), &
        (dccx(n,i),i=1,2),(dcfx(n,i),i=1,3)
      enddo
    endif
  else
    do n=1,nem
      write(it,"(//,3x,'element no. =', i3,/)") n
      if(ntype /= 0)then
        write(it,450) pr(n),se(n),sl(n),sa(n),si(n),cs(n),sn(n), &
        hf(n),vf(n),pf(n),xb(n),cnt(n),snt(n),(nod(n,i),i=1,2)
      else
        write(it,470) se(n),sl(n),sa(n),cs(n),sn(n),hf(n),(nod(n,i),i=1,2)
      endif
    enddo
  endif
! _______________________________________________
! |                                                |
! |           p r o c e s s o r   u n i t          |
! |_______________________________________________|

! time marching scheme begins here. for item=2, initial conditions
! on second derivatives of the solution are computed in the program

  if(item /= 0)then
    if(item == 1)then
      nt=nt+1
      time=time+dt
    endif
  endif

  if(item >= 3) nhbw=neq

!        initialize global matrices and vectors

  do
    do i=1,neq
      glf(i)=0.0
      do j=1,nhbw
        if(item >= 3)then
          glm(i,j)=0.0
        endif
        glk(i,j)=0.0
      enddo
    enddo

  ! do-loop for element calculations and assembly

    do ne = 1, nem
      if(model /= 4)then
        if(icont /= 1) then
          ax0=dcax(ne,1)
          ax1=dcax(ne,2)
          bx0=dcbx(ne,1)
          bx1=dcbx(ne,2)
          cx0=dccx(ne,1)
          cx1=dccx(ne,2)
          fx0=dcfx(ne,1)
          fx1=dcfx(ne,2)
          fx2=dcfx(ne,3)
        endif
        l=0
        do i=1,npe
          ni=nod(ne,i)
          if(icont == 1)then
            elx(i)=glx(ni)
          else
            elx(1)=0.0
            elx(2)=0.5*glx(ne)
            elx(npe)=glx(ne)
          endif
          if(item == 1 .or. item == 2)then
            li=(ni-1)*ndf
            do j=1,ndf
              li=li+1
              l=l+1
              elu(l)=gu0(li)
              if(item == 2 .and. nt > 0)then
                elv(l)=gu1(li)
                ela(l)=gu2(li)
              endif
            enddo
          endif
        enddo
        call coeffcnt(ielem,item,model,ndf,npe,time,ntype,ne,f3,mxelm)
      else
        call transfrm(mxelm,ne,ntype,pr,se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
      endif

      if(nprnt  /= 0)then
        nn = npe*ndf
        if(nprnt  <= 2)then
          if(ne <= 5 .and. nt <= 1)then
            write(it,"(/,3x,'element coefficient matrix, [elk]:',/)")
            do  i=1,nn
              write(it,"(2x,5e13.5)") (elk(i,j),j=1,nn)
            enddo
            if(item >= 3)then
              write(it,"(/,3x,'element coefficient matrix, [elm]:',/)")
              do i=1,nn
                write(it,"(2x,5e13.5)") (elm(i,j),j=1,nn)
              enddo
            else
              write(it,"(/,3x,'element source vector, {elf}:',/)")
              write(it,"(2x,5e13.5)") (elf(i),i=1,nn)
            endif
          endif
        endif
      endif

      !   assemble element matrices
      call assemble(nod,mxelm,mxneq,ndf,npe,ne,item,glk,glm,glf)
    enddo

    ! call subroutine constrnt to impose constraint boundary conditions,
    ! for example, inclined support conditions

    if(model == 4)then
      if(ncon /= 0)then
        call constrnt(neq,nhbw,ndf,ncon,icon,vcon,glk,glm,glf,trm,mxneq)
      endif
    endif
!!! continuar desde  aqui
    ! impose multi-point constraints using the penalty method

    if(nmpc /= 0)then
      if(nprnt == 2)then
        write(it,"(/,3x,'global coefficient matrix, [glk]:',/)")
        do i=1,neq
          write(it,"(2x,5e13.5)") (glk(i,j),j=1,nhbw)
        enddo
      endif
      vmax=0.0
      do i=1,neq
        do j=i,nhbw
          value1=dabs(glk(i,j))
          if(value1 > vmax)then
            vmax=value1
          endif
        enddo
      enddo
      pnlty=vmax*1.0e4
      do nc=1,nmpc
        ndof1=(imc1(nc,1)-1)*ndf+imc1(nc,2)
        ndof2=(imc2(nc,1)-1)*ndf+imc2(nc,2)
        glk(ndof1,1)=glk(ndof1,1)+pnlty*vmpc(nc,1)*vmpc(nc,1)
        glk(ndof2,1)=glk(ndof2,1)+pnlty*vmpc(nc,2)*vmpc(nc,2)
        glf(ndof1)=glf(ndof1)+pnlty*vmpc(nc,1)*vmpc(nc,3)
        glf(ndof2)=glf(ndof2)+pnlty*vmpc(nc,2)*vmpc(nc,3)
        if(ndof1 > ndof2)then
          nw=ndof1-ndof2+1
          glk(ndof2,nw)=glk(ndof2,nw)+pnlty*vmpc(nc,1)*vmpc(nc,2)
          glf(ndof1)=vmpc(nc,4)
        else
          nw=ndof2-ndof1+1
          glk(ndof1,nw)=glk(ndof1,nw)+pnlty*vmpc(nc,1)*vmpc(nc,2)
          glf(ndof2)=vmpc(nc,4)
        endif
      enddo
    endif

    if(nprnt == 2)then
    ! print assembled coefficient matrices if required
      write(it,"(/,3x,'global coefficient matrix, [glk]:',/)")
      do i=1,neq
        write(it,"(2x,5e13.5)") (glk(i,j),j=1,nhbw)
      enddo
      if(item >= 3)then
        write(it,"(/,3x,'global coefficient matrix, [glm]:',/)")
        do i=1,neq
          write(it,"(2x,5e13.5)") (glm(i,j),j=1,nhbw)
        enddo
      else
        write(it,"(/,3x,'global source vector, {glf}:',/)")
        write(it,"(2x,5e13.5)") (glf(i),i=1,neq)
      endif
    endif

    ! call subroutine boundary to impose essential, natural and newton's
    ! type boundary conditions on the primary and secondary variables.

    call boundary(neq,neqr,nhbw,nspv,nssv,nnbc,ndf,dt,item,alfa,ibdy,&
                  ispv,issv,inbc,uref,vspv,vssv,vnbc,glk,glm,glf,gu0,&
                  mxebc,mxnbc,mxmbc,mxneq)

    if(nprnt == 2)then
      ! print assembled coefficient matrices if required
      write(it,"(/,3x,'global coefficient matrix, [glk]:',/)")
      do i=1,neq
        write(it,"(2x,5e13.5)") (glk(i,j),j=1,nhbw)
      enddo
    endif

    if(item >= 3)then
      ! call egnsolvr to solve for the eigenvalues and eigenvectors
      call egnsolvr(neqr,glk,glm,egnval,egnvec,jvec,nrot,mxneq)
      write(it,"(/,5x,'number of rotations taken in jacobi =',i2,/)") nrot
      do nvec=1,neqr
        frqncy=dsqrt(egnval(nvec))
        write(it,"(/,5x,'eigenvalue(',i2,') = ',e14.6,2x,'sqrt(egnval) = ' ,e13.5,/,5x,'eigenvector:')") &
                nvec,egnval(nvec),frqncy
        write(it,"(2x,5e13.5)")(egnvec(i,nvec),i=1,neqr)
      enddo
      stop
    endif

    ires = 0

    ! call subroutine eqnsolvr to solve the finite-element equations

    call eqnsolvr(mxneq,mxneq,neq,nhbw,glk,glf,ires)

    if(item == 0)then
      write(it,"(/,1x,'solution (values of pvs) at the nodes: ',/)")
      write(it,"(2x,5e13.5)") (glf(ni),ni=1,neq)
    else
      if(nt == 0)then
        do i=1,neq
          gu2(i)=glf(i)
        enddo
        nt=nt+1
        time=time+dt
        cycle
      endif

      ! compute and print current values of gu0, gu1, and gu2

      do i=1,neq
        if(item == 2)then
          acclrn=a3*(glf(i)-gu0(i))-a4*gu1(i)-a5*gu2(i)
          gu1(i)=gu1(i)+a2*gu2(i)+a1*acclrn
          gu2(i)=acclrn
          gpu(i)=gu0(i)
        else
          gpu(i)=gu0(i)
        endif
        gu0(i)=glf(i)
      enddo

      diff=0.0
      soln=0.0
      do i=1,neq
        soln=soln+gu0(i)*gu0(i)
        diff=diff+(glf(i)-gpu(i))**2
      enddo

      prcnt=dsqrt(diff/soln)
      if(prcnt <= tolrns)then
        write(it,640)
        write(it,"(2x,5e13.5)") (gpu(i),i=1,neq)
        write(it,"(2x,5e13.5)") (gu0(i),i=1,neq)
        stop
      else
        if(intvl <= 0)intvl=1
        nten=(nt/intvl)*intvl
        if(nten == nt)then
          write(it,"(/,1x,'time =',e12.4,5x,'time step number =',i3,/)") time, nt
          write(it,"(/,1x,'solution (values of pvs) at the nodes: ',/)")
          write(it,"(2x,5e13.5)") (gu0(i),i=1,neq)
          if(item /= 1) then
            if(nprnt < 4)then
              write(it,"(/,2x,'first time derivative of the primary variables:',/)")
              write(it,"(2x,5e13.5)") (gu1(i),i=1,neq)
              write(it,"(/,2x,'second time derivative of the primary variables:',/)")
              write(it,"(2x,5e13.5)") (gu2(i),i=1,neq)
            endif
          endif
          nt=nt+1
          time=time+dt
        else
          nt=nt+1
          time=time+dt
          cycle
        endif
      endif
    endif

    if(nmpc == 0)then
      if(nprnt <= 1)then
        if(model == 1)then
          write(it,"(2x,55('_'),/)")
        else
          if(model == 4)then
            write(it,630)
          endif
          write(it,"(2x,78('_'),/)")
        endif

        if(model == 1)then
          write(it,"(3x,'x is the global coord. if icont=1 and it is the local', ' coord. if icont=0')")
          if(ntype == 0)then
            write(it,"(7x,' x ',5x, 'p. variable',2x,'s. variable')")
          else
            write(it,"(7x,' x ',5x, 'displacement',2x,'radial stress',2x,'hoop stress')")
          endif
        endif

        if(model == 2 .or. model == 3)then
          write(it,"(3x,'x is the global coord. if icont=1 and it is the local', ' coord. if icont=0')")
          if(ntype == 0)then
            write(it,"(7x,' x ',6x, 'deflect.',5x,'rotation',5x,'b. moment', 3x,'shear force')")
          else
            write(it,660)
          endif
        endif

        if(model == 4)then
          if(ntype == 0)then
            write(it,"(3x, 'ele force, h1      force, v1   force, h2 force, v2')")
          else
            write(it,"(3x, 'ele force, h1      force, v1 moment, m1 force, h2 force, v2 moment, m2')")
          endif
        endif

        if(model == 1)then
          write(it,"(2x,55('_'),/)")
        else
          write(it,"(2x,78('_'),/)")
        endif

        if(model <= 3)then
          call postproc(dcax,dcbx,dccx,f3,glf,glx,nod,icont,ielem,npe,&
          model,ntype,item,mxelm,mxneq,mxnod,nem,ndf)
        else
          call reaction(mxelm,mxneq,ndf,nem,nod,npe,ntype,pr,glf,&
          se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
        endif

        if(model == 1)then
          write(it,"(2x,55('_'),/)")
        else
          write(it,"(2x,78('_'),/)")
        endif
      endif
    else
      ! calculate the reactions at the points where constraints are imposed
      do nc=1,nmpc
        ndof1=(imc1(nc,1)-1)*ndf+imc1(nc,2)
        ndof2=(imc2(nc,1)-1)*ndf+imc2(nc,2)
        gu0(nc)=-pnlty*vmpc(nc,1)*(vmpc(nc,1)*glf(ndof1)+vmpc(nc,2)*glf(ndof2)-vmpc(nc,3))
        gu1(nc)=-pnlty*vmpc(nc,2)*(vmpc(nc,1)*glf(ndof1)+vmpc(nc,2)*glf(ndof2)-vmpc(nc,3))
      enddo

      write(it,"(/,3x,'forces at the constrained points:',/)")
      write(it,"(2x,5e13.5)")(gu0(i),i=1,nmpc)
      write(it,"(2x,5e13.5)")(gu1(i),i=1,nmpc)
    endif

    if(item == 0) stop
    if(nt < ntime)then
      if(prcnt > tolrns)then
        cycle
      endif
    else
      write(it,"(/,5x,'***** number of time steps exceeded ntime *****',/)")
      exit
    endif
  enddo
  close(in)
  close(it)

    !!!!!!!!!!!!!!! Deallocate dnamically allocated variables !!!!!!!!!!!!!!
    deallocate (ispv)
    deallocate (ibdy)
    deallocate (issv)
    deallocate (inbc)
    deallocate (imc1)
    deallocate (imc2)
    deallocate (nod)
    deallocate (icon)

    deallocate ( dcax )
    deallocate ( dcbx )
    deallocate ( dccx )
    deallocate ( dcfx )

    deallocate ( gu0 )
    deallocate ( gu1 )
    deallocate ( gu2 )
    deallocate ( gpu )
    deallocate ( dx  )

    deallocate ( vcon )
    deallocate ( glf )
    deallocate ( glx )

    deallocate ( glm )
    deallocate ( trm )
    deallocate ( vmpc )

    deallocate ( cs )
    deallocate ( sn )
    deallocate ( cnt)
    deallocate ( snt )
    deallocate ( xb )

    deallocate ( egnvec )
    deallocate ( glk )

    deallocate ( hf )
    deallocate ( vf )
    deallocate ( pf )
    deallocate ( f3 )

    deallocate ( pr )
    deallocate ( se )
    deallocate ( sl )
    deallocate ( sa )
    deallocate ( si )
    deallocate ( egnval )

    deallocate ( uref )
    deallocate ( vspv )
    deallocate ( vssv )
    deallocate ( vnbc )


    !!!!!!!!!!!!!!! Deallocate dnamically allocated variables !!!!!!!!!!!!!!

  stop
end ! end of fem1d Program !
