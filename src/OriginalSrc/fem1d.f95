!       program name: fem1d           length(includinmg blanks):2440 lines
!
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       *                        program fem1d                          *
!       *          (a finite element analysis computer program)         *
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         _______________________________________________________________
!       |                                                                |
!       | this is a finite element computer program for the analysis     |
!       | of the following three model equations and others:             |
!       |                                                                |
!       | 1. heat transfer, fluid mechanics, bars, and cables:           |
!       |                                                                |
!       |             ct.u* + ct.u** - (ax.u')' + cx.u = fx              |
!       |                                                                |
!       | 2. the timoshenko beam and circular plate theory:              |
!       |                                                                |
!       |              ct0.w** - [ax.(w' + s)]' + cx.w = fx              |
!       |              ct1.s** - (bx.s')' + ax.(w' + s) = 0              |
!       |                                                                |
!       | 3. the euler-bernoulli beam and circular plate theory:         |
!       |                                                                |
!       |                ct.w** + (bx.w'')'' + cx.w = fx                 |
!       |                                                                |
!       | in the above equations (') and (*) denote differentiations     |
!       | with respect to space x and time t, and ax, bx, cx, ct, and    |
!       | fx are functions of x only:                                    |
!       |                                                                |
!       |     ax = ax0 + ax1.x, bx = bx0 + bx1.x, cx = cx0 + cx1.x       |
!       |          ct = ct0 + ct1.x, fx = fx0 + fx1.x + fx2.x.x          |
!       |                                                                |
!       |     in addition to the three model equations, other equations  |
!       | (for example, disks, trusses, and frames) can be analyzed by   |
!       | the program.                                                   |
!       |________________________________________________________________|
!
!        _________________________________________________________________
!        .                                                               .
!        .               key variables used in the program               .
!        . see table 7.3.2 of the book for a description of the variables.
!        .                                                               .
!        . ndf....... number of degrees of freedom per node              .
!        . neq....... number of equations in the model (before b. c.)    .
!        . ngp....... number of gauss points used in the evaluation of   .
!        .            the element coefficients, elk , elf , elm          .
!        . nhbw...... half bandwidth of global coefficient matrix glk    .
!        . nn ....... number of total degrees of freedom in the element .
!        . npe....... number of nodes per element                        .
!        _________________________________________________________________
!        _________________________________________________________________
!        .                                                               .
!        .          dimensions of various arrays in the program          .
!        .                                                               .
!        . values of mxelm,mxnod, etc. in the parameter statement should .
!        .         be changed to meet the requirements of the problem:   .
!        .                                                               .
!        . mxelm..... maximum number of elements in the mesh:            .
!        . mxebc..... maximum number of speci. primary deg. of freedom   .
!        . mxmbc..... maximum number of speci. mixed boundary conditions .
!        . mxnbc..... maximum number of speci. secondary deg. of freedom .
!        . mxneq..... maximum number of equations in the fe model        .
!        . mxnod..... maximum number of nodes in the mesh                .
!        .                                                               .
!        . note: the following dimension statement in subroutine jacobi .
!        .        should be modified when mxneq is greater than 500:     .
!        .         dimension v(500,500),vt(500,500),w(500,500),ih(500) .
!        .        the value of mxneq should be used in place of ‘500'    .
!        _________________________________________________________________
!        .                                                               .
!        .                subroutines used in the program                .
!        .                                                               .
!        . assemble, boundary, coeffcnt, constrnt, echodata, eqnsolvr, .
!        .    eignslvr, jacobi, matrxmlt, mesh1d, postproc, reaction,    .
!        .             shape1d, timforce, timstres, transfrm             .
!        _________________________________________________________________
!
program fem1d
Use, intrinsic :: iso_fortran_env, Only : iostat_end
implicit real*8(a-h,o-z)
!implicit none
  integer           :: argc
  character(len=32) :: argv

            parameter (mxelm=250,mxneq=500,mxebc=20,mxnbc=20,mxmbc=20,&
                              mxnod=250,mxmpc=5)
            dimension dcax(mxelm,2),dcbx(mxelm,2),dccx(mxelm,2),dcfx(mxelm,3)
            dimension gu0(mxneq),gu1(mxneq),gu2(mxneq),gpu(mxneq),dx(mxnod)
            dimension ibdy(mxebc),ispv(mxebc,2),issv(mxnbc,2),inbc(mxmbc,2)
            dimension imc1(mxmpc,2),imc2(mxmpc,2),vmpc(mxmpc,4)
            dimension icon(9),vcon(9),trm(mxneq,mxneq)
            dimension glm(mxneq,mxneq),glf(mxneq),glx(mxnod),nod(mxelm,4)
            dimension cs(mxelm),sn(mxelm),cnt(mxelm),snt(mxelm),xb(mxelm)
          dimension egnval(mxneq),egnvec(mxneq,mxneq),glk(mxneq,mxneq)
          dimension pr(mxelm),se(mxelm),sl(mxelm),sa(mxelm),si(mxelm)
          dimension hf(mxelm),vf(mxelm),pf(mxelm),f3(mxelm),title(20)
          dimension uref(mxmbc),vspv(mxebc),vssv(mxnbc),vnbc(mxmbc)
          common/stf1/elk(9,9),elm(9,9),elf(9),elx(4),elu(9),elv(9),ela(9)
          common/stf2/a1,a2,a3,a4,a5,ax0,ax1,bx0,bx1,cx0,cx1,ct0,ct1,fx0,fx1,fx2
          common/io/in,it
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

         read(in,300) title
         read(in,*) model,ntype,item
         read(in,*) ielem,nem
         read(in,*) icont,nprnt

         if(model.ge.3)then
            npe=2
            if(model.eq.4 .and. ntype.ge.1)then
               ndf=3
            else
               ndf=2
            endif
            if(model.eq.4 .and. ntype.eq.2)then
               ielem=1
            else
               ielem=0
            endif
         else
            if(model.eq.2)then
               ndf=2
               if(ntype.gt.1)ielem=1
            else
               ndf=1
            endif
            npe=ielem+1
         endif

!       data input for bar-like and beam problems (model=1,2, and 3)

           if(model.ne.4)then
              if(icont.ne.0)then
                 nnm = nem*(npe-1)+1
                 nem1=nem + 1
                 read(in,*) (dx(i), i=1,nem1)
                 call mesh1d(nem,npe,nod,mxelm,mxnod,dx,glx)
                 read(in,*) ax0,ax1
                 read(in,*) bx0,bx1
                 read(in,*) cx0,cx1
                 if(item.lt.3)then
                    read(in,*) fx0,fx1,fx2
                 endif
              else

!       read glx, nod, and element-wise continuous coefficients [dc.x]

                 read(in,*)nnm
                 do 10 n=1,nem
                 read(in,*) (nod(n,i),i=1,npe), glx(n)
                 read(in,*) (dcax(n,i),i=1,2)
                 read(in,*) (dcbx(n,i),i=1,2)
                 read(in,*) (dccx(n,i),i=1,2)
        10       read(in,*) (dcfx(n,i),i=1,3)
              endif
           else

!       input data for plane truss or frame structures (model=4)

               read(in,*)nnm
               if(ntype.ne.0)then
                  do 20 n=1,nem
                  read(in,*) pr(n),se(n),sl(n),sa(n),si(n),cs(n),sn(n)
                  read(in,*) hf(n),vf(n),pf(n),xb(n),cnt(n),snt(n)
        20        read(in,*) (nod(n,i),i=1,2)
               else
                  do 30 n=1,nem
                  read(in,*) se(n),sl(n),sa(n),cs(n),sn(n),hf(n)
        30        read(in,*) (nod(n,i),i=1,2)
               endif
                                   read(in,*) ncon
               if(ncon.ne.0)then
                  do 35 i=1, ncon
        35        read(in,*) icon(i),vcon(i)
               endif
           endif
           neq=nnm*ndf

!       read data on boundary conditions of three kinds: dirichlet (pv)
!                 neumann (sv), and newton's (mixed) types

       read(in,*) nspv
       if(nspv.ne.0)then
          do 40 nb=1,nspv
          if(item.gt.2)then
             read(in,*) (ispv(nb,j),j=1,2)
          else
             read(in,*) (ispv(nb,j),j=1,2),vspv(nb)
          endif
    40    continue
       endif

       if(item.le.2)then
          read(in,*) nssv
          if(nssv.ne.0)then
             do 50 ib=1,nssv
    50       read(in,*) (issv(ib,j),j=1,2),vssv(ib)
          endif
       endif

       read(in,*) nnbc
       if(nnbc.ne.0)then
          do 60 i=1, nnbc
    60    read(in,*) (inbc(i,j),j=1,2),vnbc(i),uref(i)
       endif

!       read data on multi-point constraints

       read(in,*) nmpc
       if(nmpc.ne.0)then
          do 65 i=1, nmpc
    65    read(in,*)(imc1(i,j),j=1,2),(imc2(i,j),j=1,2),(vmpc(i,j),j=1,4)
       endif

        if(item .ne. 0)then

!       input data here for time-dependent problems

            if(item.le.3)then
               read(in,*) ct0,ct1
            endif
            if(item.le.2)then
               read(in,*) dt,alfa,gama
               read(in,*) incond,ntime,intvl
               a1=alfa*dt
               a2=(1.0-alfa)*dt
               if(incond.ne.0)then
                  read(in,*) (gu0(i),i=1,neq)
                  else
                     do 70 i=1,neq
        70           gu0(i)=0.0
                  endif
                  if(item.eq.2)then
                     a3=2.0/gama/(dt*dt)
                     a4=a3*dt
                     a5=1.0/gama-1.0
                     if(incond.ne.0)then
                        read(in,*) (gu1(i),i=1,neq)
                     else
                        do 80 i=1,neq
                        gu1(i)=0.0
        80              gu2(i)=0.0
                     endif
                  endif
              endif
           endif

!       ----------------------------------------------------------------
!        e n d       o f      t h e       i n p u t          d a t a
!       ----------------------------------------------------------------

!       compute the half bandwidth of the coefficient matrix glk

           nhbw=0.0
           do 90 n=1,nem
           do 90 i=1,npe
           do 90 j=1,npe
           nw=(iabs(nod(n,i)-nod(n,j))+1)*ndf
        90 if(nhbw.lt.nw) nhbw=nw

!       ----------------------------------------------------------------
!                 p r i n t    t h e    i n p u t    d a t a
!       ----------------------------------------------------------------

           write(it,530)
           write(it,310)
           write(it,530)
           write(it,300) title
           write(it,320) model,ntype
           write(it,350) ielem,ndf,nem,neq,nhbw,nspv,nssv,nnbc,nmpc

           if(item.ne.0)then
              if(item.le.2)then
                 write(it,330)
                 write(it,390) ct0,ct1,alfa,gama,dt,ntime,intvl
                 if(incond.ne.0)then
                    write(it,370)
                    write(it,540) (gu0(i),i=1,neq)
                 if(item.eq.2)then
                    write(it,380)
                    write(it,540) (gu1(i),i=1,neq)
                 endif
              endif
           else
              write(it,340)
              if(item.le.3)then
                 write(it,400) ct0,ct1
              endif
           endif
        endif

        if(nspv.ne.0)then
           write(it,480)
           do 100 ib=1,nspv
           if(item.le.2)then
              write(it,490) (ispv(ib,j),j=1,2),vspv(ib)
           else
              write(it,490) (ispv(ib,j),j=1,2)
           endif
    100    continue
        endif

        if(nssv.ne.0)then
           write(it,500)
           do 110 ib=1,nssv
    110    write(it,490) (issv(ib,j),j=1,2),vssv(ib)
        endif

        if(nnbc.ne.0)then
           write(it,510)
           do 120 i=1,nnbc
    120    write(it,490) (inbc(i,j),j=1,2),vnbc(i),uref(i)
        endif

         if(nmpc.ne.0)then
            write(it,515)
            do 125 i=1, nmpc
    125     write(it,495)(imc1(i,j),j=1,2),(imc2(i,j),j=1,2),(vmpc(i,j),j=1,4)
         endif

        if(model.ne.4)then
           if(icont.eq.1)then
              write(it,410)
              write(it,540) (glx(i),i=1,nnm)
              write(it,420)
              if(model.ne.3)then
                 write(it,440) ax0,ax1,bx0,bx1,cx0,cx1,fx0,fx1,fx2
                    else
                         write(it,445) ax0,ax1,bx0,bx1,cx0,cx1
                    endif
                 else
                    do 130 n=1,nem
                    write(it,430) n,glx(n)
        130         write(it,440) (dcax(n,i),i=1,2),(dcbx(n,i),i=1,2),(dccx(n,i),i=1,2),(dcfx(n,i),i=1,3)

                 endif
              else
                 do 140 n=1,nem
                 write(it,460) n
                 if(ntype.ne.0)then
                    write(it,450) pr(n),se(n),sl(n),sa(n),si(n),cs(n),sn(n), &
                                    hf(n),vf(n),pf(n),xb(n),cnt(n),snt(n), &
                                    (nod(n,i),i=1,2)
                 else
                    write(it,470) se(n),sl(n),sa(n),cs(n),sn(n),hf(n), &
                                    (nod(n,i),i=1,2)
                 endif
        140      continue
              endif
!                    _______________________________________________
!                  |                                                |
!                  |           p r o c e s s o r   u n i t          |
!                  |_______________________________________________|

!        time marching scheme begins here. for item=2, initial conditions
!        on second derivatives of the solution are computed in the program

            if(item.ne.0)then
               if(item.eq.1)then
                  nt=nt+1
                  time=time+dt
               endif
            endif

            if(item.ge.3)nhbw=neq

!        initialize global matrices and vectors

        150 do 160 i=1,neq
            glf(i)=0.0
            do 160 j=1,nhbw
            if(item.ge.3)then
               glm(i,j)=0.0
            endif
        160 glk(i,j)=0.0

!         do-loop for element calculations and assembly

          do 200 ne = 1, nem
          if(model.ne.4)then
             if(icont.ne.1) then
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
             do 180 i=1,npe
             ni=nod(ne,i)
             if(icont.eq.1)then
                 elx(i)=glx(ni)
             else
                 elx(1)=0.0
                 elx(2)=0.5*glx(ne)
                 elx(npe)=glx(ne)
             endif
             if(item.eq.1 .or. item.eq.2)then
                 li=(ni-1)*ndf
                 do 170 j=1,ndf
                 li=li+1
                 l=l+1
                 elu(l)=gu0(li)
                 if(item.eq.2 .and. nt.gt.0)then
                    elv(l)=gu1(li)
                    ela(l)=gu2(li)
                 endif
    170          continue
             endif
    180      continue

             call coeffcnt(ielem,item,model,ndf,npe,time,ntype,ne,f3,mxelm)
          else
             call transfrm(mxelm,ne,ntype,pr,se,sl,sa,si,cs,sn,cnt,snt,&
                           hf,vf,pf,xb)
          endif

            if(nprnt .ne.0)then
               nn = npe*ndf
               if(nprnt .le.2)then
                  if(ne.le.5 .and. nt.le.1)then
                     write(it,550)
                     do 190 i=1,nn
        190          write(it,540) (elk(i,j),j=1,nn)
                     if(item.ge.3)then
                        write(it,360)
                        do 195 i=1,nn
        195             write(it,540) (elm(i,j),j=1,nn)
                     else
                        write(it,560)
                        write(it,540) (elf(i),i=1,nn)
                     endif
                  endif
               endif
            endif

!        assemble element matrices

            call assemble(nod,mxelm,mxneq,ndf,npe,ne,item,glk,glm,glf)

        200 continue

!        call subroutine constrnt to impose constraint boundary conditions,
!        for example, inclined support conditions

            if(model.eq.4)then
               if(ncon.ne.0)then
            call constrnt(neq,nhbw,ndf,ncon,icon,vcon,glk,glm,glf,trm,mxneq)
               endif
            endif

!        impose multi-point constraints using the penalty method

            if(nmpc.ne.0)then
               if(nprnt.eq.2)then
                  write(it,570)
                  do 201 i=1,neq
        201       write(it,540) (glk(i,j),j=1,nhbw)
               endif
               vmax=0.0
               do 204 i=1,neq
               do 204 j=i,nhbw
               value=dabs(glk(i,j))
               if(value.gt.vmax)then
                  vmax=value
               endif
    204    continue
           pnlty=vmax*1.0e4
           do 205 nc=1,nmpc
           ndof1=(imc1(nc,1)-1)*ndf+imc1(nc,2)
           ndof2=(imc2(nc,1)-1)*ndf+imc2(nc,2)
           glk(ndof1,1)=glk(ndof1,1)+pnlty*vmpc(nc,1)*vmpc(nc,1)
           glk(ndof2,1)=glk(ndof2,1)+pnlty*vmpc(nc,2)*vmpc(nc,2)
           glf(ndof1)=glf(ndof1)+pnlty*vmpc(nc,1)*vmpc(nc,3)
           glf(ndof2)=glf(ndof2)+pnlty*vmpc(nc,2)*vmpc(nc,3)
           if(ndof1.gt.ndof2)then
              nw=ndof1-ndof2+1
              glk(ndof2,nw)=glk(ndof2,nw)+pnlty*vmpc(nc,1)*vmpc(nc,2)
              glf(ndof1)=vmpc(nc,4)
           else
              nw=ndof2-ndof1+1
              glk(ndof1,nw)=glk(ndof1,nw)+pnlty*vmpc(nc,1)*vmpc(nc,2)
              glf(ndof2)=vmpc(nc,4)
           endif
    205    continue
        endif

            if(nprnt.eq.2)then

!           print assembled coefficient matrices if required


           write(it,570)
           do 210 i=1,neq
    210    write(it,540) (glk(i,j),j=1,nhbw)
           if(item.ge.3)then
              write(it,575)
              do 215 i=1,neq
    215       write(it,540) (glm(i,j),j=1,nhbw)
           else
              write(it,580)
              write(it,540) (glf(i),i=1,neq)
           endif
        endif

!           call subroutine boundary to impose essential, natural and newton's
!           type boundary conditions on the primary and secondary variables.

            call boundary(neq,neqr,nhbw,nspv,nssv,nnbc,ndf,dt,item,alfa,ibdy,&
                          ispv,issv,inbc,uref,vspv,vssv,vnbc,glk,glm,glf,gu0,&
                          mxebc,mxnbc,mxmbc,mxneq)

            if(nprnt.eq.2)then

!           print assembled coefficient matrices if required

               write(it,570)
               do 211 i=1,neq
        211    write(it,540) (glk(i,j),j=1,nhbw)
            endif

              if(item.ge.3)then

!          call egnsolvr to solve for the eigenvalues and eigenvectors

                 call egnsolvr(neqr,glk,glm,egnval,egnvec,jvec,nrot,mxneq)

               write(it,690) nrot
               do 230 nvec=1,neqr
               frqncy=dsqrt(egnval(nvec))
               write(it,700)nvec,egnval(nvec),frqncy
        230    write(it,540)(egnvec(i,nvec),i=1,neqr)
               stop
            endif

              ires = 0

!          call subroutine eqnsolvr to solve the finite-element equations

              call eqnsolvr(mxneq,mxneq,neq,nhbw,glk,glf,ires)

            if(item.eq.0)then
               write(it,590)
               write(it,540) (glf(ni),ni=1,neq)
            else
               if(nt.eq.0)then
                  do 240 i=1,neq
        240       gu2(i)=glf(i)
                  nt=nt+1
                  time=time+dt
                  goto 150
               endif

!          compute and print current values of gu0, gu1, and gu2

                 do 250 i=1,neq
                 if(item.eq.2)then
                    acclrn=a3*(glf(i)-gu0(i))-a4*gu1(i)-a5*gu2(i)
                    gu1(i)=gu1(i)+a2*gu2(i)+a1*acclrn
                    gu2(i)=acclrn
                    gpu(i)=gu0(i)
                 else
                    gpu(i)=gu0(i)
                 endif
        250      gu0(i)=glf(i)

            diff=0.0
            soln=0.0
            do 260 i=1,neq
            soln=soln+gu0(i)*gu0(i)
    260     diff=diff+(glf(i)-gpu(i))**2
            prcnt=dsqrt(diff/soln)
            if(prcnt.le.tolrns)then
               write(it,640)
               write(it,540) (gpu(i),i=1,neq)
               write(it,540) (gu0(i),i=1,neq)
               stop
            else
               if(intvl.le.0)intvl=1
                  nten=(nt/intvl)*intvl
                  if(nten.eq.nt)then
                     write(it,600) time, nt
                     write(it,590)
                     write(it,540) (gu0(i),i=1,neq)
                     if(item.ne.1) then
                        if(nprnt.lt.4)then
                           write(it,645)
                           write(it,540) (gu1(i),i=1,neq)
                           write(it,646)
                           write(it,540) (gu2(i),i=1,neq)
                        endif
                     endif
                     nt=nt+1
                     time=time+dt
                  else
                     nt=nt+1
                     time=time+dt
                     goto 150
                  endif
               endif
            endif

          if(nmpc.eq.0)then
             if(nprnt.le.1)then
                if(model.eq.1)then
                   write(it,530)
                else
                   if(model.eq.4)then
                      write(it,630)
                   endif
                   write(it,520)
                endif

                  if(model.eq.1)then
                     write(it,647)
                     if(ntype.eq.0)then
                        write(it,610)
                     else
                        write(it,620)
                     endif
                  endif

                  if(model.eq.2 .or. model.eq.3)then
                     write(it,647)
                     if(ntype.eq.0)then
                        write(it,650)
                     else
                        write(it,660)
                     endif
                  endif

                  if(model.eq.4)then
                     if(ntype.eq.0)then
                        write(it,680)
                     else
                        write(it,670)
                     endif
                  endif

                  if(model.eq.1)then
                     write(it,530)
                  else
                     write(it,520)
                  endif

                  if(model.le.3)then
                  call postproc(dcax,dcbx,dccx,f3,glf,glx,nod,icont,ielem,npe,&
                                model,ntype,item,mxelm,mxneq,mxnod,nem,ndf)
                  else
                  call reaction(mxelm,mxneq,ndf,nem,nod,npe,ntype,pr,glf,&
                                se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
                  endif

               if(model.eq.1)then
                  write(it,530)
               else
                  write(it,520)
               endif
            endif
         else

!        calculate the reactions at the points where constraints are imposed

               do 280 nc=1,nmpc
               ndof1=(imc1(nc,1)-1)*ndf+imc1(nc,2)
               ndof2=(imc2(nc,1)-1)*ndf+imc2(nc,2)
           gu0(nc)=-pnlty*vmpc(nc,1)*(vmpc(nc,1)*glf(ndof1) &
                    +vmpc(nc,2)*glf(ndof2)-vmpc(nc,3))
           gu1(nc)=-pnlty*vmpc(nc,2)*(vmpc(nc,1)*glf(ndof1) &
                    +vmpc(nc,2)*glf(ndof2)-vmpc(nc,3))
   280     continue
           write(it,545)
           write(it,540)(gu0(i),i=1,nmpc)
           write(it,540)(gu1(i),i=1,nmpc)
        endif

       if(item.eq.0)stop
       if(nt.lt.ntime)then
          if(prcnt.gt.tolrns)then
             goto 150
          endif
       else
          write(it,710)
       endif
!      ----------------------------------------------------------------
!                        f   o   r   m   a    t  s
!      ----------------------------------------------------------------
  300 format(20a4)
  310 format(8x,'output from program    fem1d   by j n reddy')
  320 format(/,4x,'*** analysis of model',i2,', and type',i2, &
              ' problem ***',/,15x,'(see the code below)',/, &
              /,4x,'model=1,ntype=0: a problem described by model eq. 1', &
              /,4x,'model=1,ntype=1: a circular disk (plane stress) ', &
              /,4x,'model=1,ntype>1: a circular disk (plane strain) ', &
              /,4x,'model=2,ntype=0: a timoshenko beam (rie) problem', &
              /,4x,'model=2,ntype=1: a timoshenko plate (rie) problem', &
              /,4x,'model=2,ntype=2: a timoshenko beam (cie) problem', &
              /,4x,'model=2,ntype>2: a timoshenko plate (cie) problem', &
              /,4x,'model=3,ntype=0: a euler-bernoulli beam problem', &
              /,4x,'model=3,ntype>0: a euler-bernoulli circular plate', &
              /,4x,'model=4,ntype=0: a plane truss problem', &
              /,4x,'model=4,ntype=1: a euler-bernoulli frame problem', &
              /,4x,'model=4,ntype=2: a timoshenko (cie) frame problem',/)
  330 format(/,4x,'time-dependent (transient) analysis ',/)
  340 format(/,4x,'e i g e n v a l u e a n a l y s i s',/)
  350 format(/,8x, 'element type (0, hermite,>0, lagrange)..=',i4,/, &
                8x, 'no. of deg. of freedom per node, ndf....=',i4,/, &
                8x, 'no. of elements in the mesh, nem........=',i4,/, &
                8x, 'no. of total dof in the model, neq......=',i4,/, &
                8x, 'half bandwidth of matrix [glk], nhbw ...=',i4,/, &
                8x, 'no. of specified primary dof, nspv......=',i4,/, &
                8x, 'no. of specified secondary dof, nssv....=',i4,/, &
                8x, 'no. of specified newton b. c.: nnbc.....=',i4,/, &
                8x, 'no. of speci. multi-pt. cond.: nmpc.....=',i4)
  360 format(/,3x,'element coefficient matrix, [elm]:',/)
  370 format(/,3x, 'initial conditions on the primary variables:',/)
        380 format(/,3x, 'initial cond. on time der. of primary variables:',/)
        390 format(/,8x,'coefficient, ct0........................=',e12.4,/, &
                     8x,'coefficient, ct1........................=',e12.4,/, &
                     8x,'parameter, alfa.........................=',e12.4,/, &
                     8x,'parameter, gama.........................=',e12.4,/, &
                     8x,'time increment, dt......................=',e12.4,/, &
                     8x,'no. of time steps, ntime................=',i4,/, &
                     8x,'time-step interval to print soln., intvl=',i4,/)
        400 format(/,8x,'coefficient, ct0........................=',e12.4,/, &
                     8x,'coefficient, ct1........................=',e12.4,/)
        410 format(/,3x,'global coordinates of the nodes, {glx}:',/)
        420 format(/,3x,'coefficients of the differential equation:',/)
        430 format(/,5x,'properties of element =',i3,//, &
                     8x,'element length, h ....... =',e12.4)
        440 format( 8x,'ax0 =',e12.4,5x,'ax1 =',e12.4,/, &
                     8x,'bx0 =',e12.4,5x,'bx1 =',e12.4,/, &
                     8x,'cx0 =',e12.4,5x,'cx1 =',e12.4,/, &
                     8x,'fx0 =',e12.4,5x,'fx1 =',e12.4,5x,'fx2 =',e12.4,/)
        445 format( 8x,'ax0 =',e12.4,5x,'ax1 =',e12.4,/, &
                     8x,'bx0 =',e12.4,5x,'bx1 =',e12.4,/, &
                     8x,'cx0 =',e12.4,5x,'cx1 =',e12.4,/)
        450 format(8x,'the poisson ratio,           pr........ =',e12.4,/, &
                   8x,'modulus of elasticity,       se........ =',e12.4,/, &
                   8x,'length of the element,       sl........ =',e12.4,/, &
                   8x,'area of cross section,       sa........ =',e12.4,/, &
                   8x,'moment of inertia,           si........ =',e12.4,/, &
                   8x,'cosine of orientation,       cn........ =',e12.4,/, &
                   8x,'sine of orientation,         sn........ =',e12.4,/, &
                   8x,'axial body force (constant), hf........ =',e12.4,/, &
                   8x,'transverse body force (cnst),vf........ =',e12.4,/, &
                   8x,'internal point force,        pf........ =',e12.4,/, &
                   8x,'location of pf from node 1, xb........ =',e12.4,/, &
                   8x,'orientation of pf: cosine,   cst....... =',e12.4,/, &
                   8x,'orientation of pf: sine,     snt....... =',e12.4,/, &
                   8x,'nodal connectivity:          nod(i,j).. =',2i6,/)

        460 format(//,3x,'element no. =', i3,/)
        470 format(8x,'modulus of elasticity,       se........ =',e12.4,/, &
                   8x,'length of the element,       sl........ =',e12.4,/, &
                   8x,'area of cross section,       sa........ =',e12.4,/, &
                   8x,'cosine of orientation,       cn........ =',e12.4,/, &
                   8x,'sine of orientation,         sn........ =',e12.4,/, &
                   8x,'axial body force (constant), hf........ =',e12.4,/, &
                   8x,'nodal connectivity:          nod(i,j).. =',2i6,/)
        480 format(/,3x, 'boundary information on primary variables:',/)
        490 format(5x,2i5,2e15.5)
        495 format(5x,2i5,2x,2i5,/,5x,4e15.5)
        500 format(/,3x, 'boundary information on secondary variables:',/)
        510 format(/,3x, 'boundary information on mixed boundary cond.:',/)
        515 format(/,3x, 'multi-point constraint information:',/)
    520 format(2x,78('_'),/)
    530 format(2x,55('_'),/)
    540 format(2x,5e13.5)
    545 format(/,3x,'forces at the constrained points:',/)
    550 format(/,3x,'element coefficient matrix, [elk]:',/)
    560 format(/,3x,'element source vector, {elf}:',/)
    570 format(/,3x,'global coefficient matrix, [glk]:',/)
    575 format(/,3x,'global coefficient matrix, [glm]:',/)
    580 format(/,3x,'global source vector, {glf}:',/)
    590 format(/,1x,'solution (values of pvs) at the nodes: ',/)
    600 format(/,1x,'time =',e12.4,5x,'time step number =',i3,/)
    610 format(7x,' x ',5x, 'p. variable',2x,'s. variable')
    620 format(7x,' x ',5x, 'displacement',2x,'radial stress',2x, &
               'hoop stress')
    630 format(/,9x,'generalized forces in the element coordinates',/, &
         5x,'(second line gives the results in the global coordinates)')
    640 format(/,3x,'*** the solution has reached a steady state ***', &
                /,3x,'solution at the two consecutive time steps follows:')
    645 format(/,2x,'first time derivative of the primary variables:',/)
    646 format(/,2x,'second time derivative of the primary variables:',/)
    647 format(3x,'x is the global coord. if icont=1 and it is the local', &
                ' coord. if icont=0')
    650 format(7x,' x ',6x, 'deflect.',5x,'rotation',5x,'b. moment', &
               3x,'shear force')
    660 format(7x,' x ',6x, 'deflect.',5x,'rotation',4x,'moment, mr', &
               3x,'moment, mt',3x,'shear force')
    670 format(3x, 'ele force, h1      force, v1 moment, m1 force, h2 &
        force, v2 moment, m2')
    680 format(3x, 'ele force, h1      force, v1   force, h2 force, v2')
    690 format(/,5x,'number of rotations taken in jacobi =',i2,/)
    700 format(/,5x,'eigenvalue(',i2,') = ',e14.6,2x,'sqrt(egnval) = ', &
                e13.5,/,5x,'eigenvector:')
    710 format(/,5x,'***** number of time steps exceeded ntime *****',/)
         close(in)
         close(it)
         stop
         end


        subroutine assemble(nod,mxelm,mxneq,ndf,npe,ne,item,glk,glm,glf)
!       __________________________________________________________________

!         the subroutine is called in main to assemble element coefficient
!         matrices (in a upper banded matrix form) and right-hand vectors

!          {elf}.... element source vector, {f}
!          {elk}.... element coefficient matrix, [k]
!          {elm}.... element coefficient matrix, [m]
!          [nod].... connectivity matrix, [b]
!       __________________________________________________________________

                implicit real*8 (a-h,o-z)
                dimension glk(mxneq,mxneq),glm(mxneq,mxneq),glf(mxneq),&
                           nod(mxelm,4)
                common/stf1/elk(9,9),elm(9,9),elf(9),elx(4),elu(9),elv(9),ela(9)
                if(item.le.2)then

!           assemble element coefficient matrix elk and source vector elf

                  do 50 i = 1, npe
                  nr = (nod(ne,i) - 1)*ndf
                  do 40 ii = 1, ndf
                  nr = nr + 1
                  l = (i-1)*ndf + ii
                  glf(nr) = glf(nr) + elf(l)
                  do 30 j = 1, npe
                  ncl = (nod(ne,j)-1)*ndf
                  do 20 jj = 1, ndf
                  m = (j-1)*ndf + jj
                  nc = ncl-nr+jj+1
                  if(nc)20,20,10
         10       glk(nr,nc) = glk(nr,nc) + elk(l,m)
         20       continue
         30       continue
         40       continue
         50       continue
               else

!           assemble element matrices into full global matrices

                  do 100 i=1,npe
                  nr=(nod(ne,i)-1)*ndf
                  do 90 ii=1,ndf
                  nr=nr+1
                  l=(i-1)*ndf+ii
                  do 80 j=1,npe
                  nc=(nod(ne,j)-1)*ndf
                  do 70 jj=1,ndf
                  m=(j-1)*ndf+jj
                  nc=nc+1
                  glk(nr,nc)=glk(nr,nc)+elk(l,m)
         60       glm(nr,nc)=glm(nr,nc)+elm(l,m)
         70       continue
         80       continue
         90       continue
        100       continue

               endif
               return
               end
         subroutine boundary(neq,neqr,nhbw,nspv,nssv,nnbc,ndf,dt,item,alfa, &
                             ibdy,ispv,issv,inbc,uref,vspv,vssv,vnbc, &
                             glk,glm,glf,gu0,mxebc,mxnbc,mxmbc,mxneq)
!        __________________________________________________________________

!           the subroutine is called in main to implement specified boundary
!            conditions on the assembled system of finite element equations
!           __________________________________________________________________

         implicit real*8 (a-h,o-z)
         dimension ispv(mxebc,2),issv(mxnbc,2),inbc(mxmbc,2),ibdy(mxebc)
         dimension uref(mxmbc),vspv(mxebc),vssv(mxnbc),vnbc(mxmbc)
         dimension glk(mxneq,mxneq),glm(mxneq,mxneq),glf(mxneq),gu0(mxneq)

!           impose boundary conditions for static and time-dependent problems

         if(item.le.2)then

!           include specified primary degrees of freedom

               if(nspv.ne.0)then
                  do 30 nb = 1,nspv
                  ie=(ispv(nb,1)-1)*ndf+ispv(nb,2)
                  it=nhbw-1
                  i=ie-nhbw
                  do 10 ii=1,it
                  i=i+1

                  if(i .ge. 1)then
                     j=ie-i+1
                     glf(i)=glf(i)-glk(i,j)*vspv(nb)
                     glk(i,j)=0.0
                  endif
    10            continue
                  glk(ie,1)=1.0
                  glf(ie)=vspv(nb)
                  i=ie
                  do 20 ii=2,nhbw
                  i=i+1
                  if(i .le. neq)then
                     glf(i)=glf(i)-glk(ie,ii)*vspv(nb)
                     glk(ie,ii)=0.0
                  endif
    20            continue
    30            continue
              endif
              if(nssv.ne.0)then

!          include specified secondary degrees of freedom

                   do 40 nf = 1,nssv
                   nb=(issv(nf,1)-1)*ndf+issv(nf,2)
                   if(item.eq.1)glf(nb)=glf(nb)+vssv(nf)*dt
         40        if(item.ne.1)glf(nb)=glf(nb)+vssv(nf)
                endif
                if(nnbc.ne.0)then

!          include specified mixed boundary conditions

                   do 50 ic=1,nnbc
                   nc=(inbc(ic,1)-1)*ndf+inbc(ic,2)
                   if(item.eq.1)then
                      glk(nc,1)=glk(nc,1)+alfa*dt*vnbc(ic)
                      glf(nc)=glf(nc)+dt*vnbc(ic)*(uref(ic) &
                                     -(1.0-alfa)*gu0(nc))
                   else
                      glk(nc,1)=glk(nc,1)+vnbc(ic)
                      glf(nc)=glf(nc)+vnbc(ic)*uref(ic)
                   endif
         50        continue
                endif
              else

!          impose boundary conditions for eigenvalue problems

                if(nnbc.ne.0)then

!          include specified mixed boundary conditions

                   do 70 ic=1,nnbc
                   nc=(inbc(ic,1)-1)*ndf+inbc(ic,2)
                   glk(nc,nc)=glk(nc,nc)+vnbc(ic)
         70        continue
                endif

!          include specified primary degrees of freedom

                if(nspv.ne.0)then
                   do 80 ib=1,nspv
         80        ibdy(ib)=(ispv(ib,1)-1)*ndf+ispv(ib,2)
                   do 120 i=1,nspv
                   imax=ibdy(i)
                   do 110 j=i,nspv
                   if(ibdy(j).ge.imax)then
                      imax=ibdy(j)
                      ikept=j
                   endif
        110        continue
                ibdy(ikept)=ibdy(i)
                ibdy(i)=imax
    120         continue
                neqr = neq
                do 180 i=1,nspv
                ib=ibdy(i)
                if(ib .lt. neqr)then
                   neqr1=neqr-1
                   do 160 ii=ib,neqr1
                   do 140 jj=1,neqr
                   glm(ii,jj)=glm(ii+1,jj)
    140            glk(ii,jj)=glk(ii+1,jj)
                   do 150 jj=1,neqr
                   glm(jj,ii)=glm(jj,ii+1)
    150            glk(jj,ii)=glk(jj,ii+1)
    160            continue
                endif
                neqr=neqr-1
    180         continue
            endif
          endif
          return
          end



          subroutine coeffcnt(ielem,item,model,ndf,npe,time,ntype,&
                              ne,f3,mxelm)
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
!            [elm]..... element 'mass' matrix [m]
!         __________________________________________________________________

            implicit real*8(a-h,o-z)
            common/stf1/elk(9,9),elm(9,9),elf(9),elx(4),elu(9),elv(9),ela(9)
            common/stf2/a1,a2,a3,a4,a5,ax0,ax1,bx0,bx1,cx0,cx1,ct0,ct1,fx0,&
                        fx1,fx2
            common/shp/sf(4),gdsf(4),gddsf(4),gj
            dimension gauspt(5,5),gauswt(5,5),f3(mxelm)

           data gauspt/5*0.0d0,-.57735027d0,.57735027d0,3*0.0d0,-.77459667d0, &
            0.0d0,.77459667d0,2*0.0d0,-.86113631d0,-.33998104d0,.33998104d0, &
           .86113631d0,0.0d0,-.906180d0,-.538469d0,0.0d0,.538469d0,.906180d0/

           data gauswt/2.0d0,4*0.0d0,2*1.0d0,3*0.0d0,.55555555d0,.88888888d0, &
            0.55555555d0,2*0.0d0,.34785485d0,2*.65214515d0,.34785485d0,0.0d0, &
            0.236927d0,.478629d0,.568889d0,.478629d0,.236927d0/

           nn=ndf*npe
           h = elx(npe) - elx(1)
           if(ielem .eq. 0)then
              ngp=4
           else
              ngp = ielem+1
           endif
           do 10 j=1,nn
           if(item.le.2)then
              elf(j) = 0.0
           endif
           do 10 i=1,nn
           if(item.gt.0)then
              elm(i,j)=0.0
           endif
        10 elk(i,j)=0.0

           if(model.ne.2)then

!       do-loop on number of gauss points begins here

              do 100 ni=1,ngp
              xi = gauspt(ni,ngp)

!       call subroutine shape1d to evaluate the interpolation functions
!            and their global derivatives at the gauss point xi

              call shape1d(h,ielem,npe,xi)
              const = gj*gauswt(ni,ngp)
              if(ielem.eq.0)then
                 x = elx(1) + 0.5*h*(1.0+xi)
              else
                 x = 0.0
                 do 30 j=1,npe
    30         x = x + sf(j)*elx(j)
            endif

!        compute coefficient matrices and vectors for vaious model problems
!            governed by single second-order and fourth-order equations
!                         (model = 1 or 3; ntype = 0 or 1)

            cx=cx0+cx1*x
            if(item.ne.3) then
               fx=fx0+fx1*x+fx2*x*x
            endif
            if(item.gt.0)then
               ct=ct0+ct1*x
            endif
            if(model.eq.1)then

!        coefficients for all single-variable problems (model=1)

                 if(ntype.eq.0)then

!        all problems governed by model equation (3.1) (ntype=0)

                    ax=ax0+ax1*x
                    do 50 j = 1,nn
                    if(item.le.2)then
                       elf(j) = elf(j) + const*sf(j)*fx
                    endif
                    do 50 i = 1,nn
                    if(item.ne.0)then
                       elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct
                    endif
                    aij = const*gdsf(i)*gdsf(j)
                    cij = const*sf(i)*sf(j)
    50              elk(i,j)=elk(i,j) + ax*aij + cx*cij
                 else

!        radially symmetric elasticity problems (model=1, ntype>0)
!               ax0=e1, ax1=e2, bx0=nu12, bx1=h, thickness

                    anu21=bx0*ax0/ax1
                    if(ntype.eq.1)then
                       c11=bx1*ax0/(1.0-bx0*anu21)
                       c22=c11*(ax1/ax0)
                       c12=bx0*c22
                    else
                       denom=1.0-bx0-anu21
                       c11=bx1*ax0*(1.0-bx0)/(1.0+bx0)/denom
                       c22=bx1*ax1*(1.0-anu21)/(1.0+anu21)/denom
                       c12=bx0*c22
                    endif
                         do 60 j=1,nn
                         if(item.le.2)then
                            elf(j) = elf(j) + const*sf(j)*fx*x
                         endif
                         do 60 i=1,nn
                         if(item.ne.0)then
                            elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct*x
                         endif
                         aij = const*gdsf(i)*gdsf(j)*c11*x
                         cij = const*sf(i)*sf(j)*cx*x
                         dij = const*(gdsf(i)*sf(j)+sf(i)*gdsf(j))*c12
                         eij = const*sf(i)*sf(j)*c22/x
        60               elk(i,j)=elk(i,j) + aij + cij + dij + eij
                      endif
               else

!         coefficients for the euler-bernoulli theory (model=2)

                      if(ntype.eq.0)then

!         the euler-bernoulli beam element (model=1 and ntype=0)

                         bx=bx0+bx1*x
                         cx=cx0+cx1*x
                         do 70 j = 1,nn
                         if(item.le.2)then
                            elf(j) = elf(j) + const*sf(j)*fx
                         endif
                         do 70 i = 1,nn
                         if(item.gt.0)then
                            if(item.le.3)then
                               elm(i,j) = elm(i,j) + const*sf(i)*sf(j)*ct
                            else
                               elm(i,j) = elm(i,j) + const*gdsf(i)*gdsf(j)
                            endif
                         endif
                         bij = const*gddsf(i)*gddsf(j)
                         cij = const*sf(i)*sf(j)
        70               elk(i,j)=elk(i,j) + bx*bij + cx*cij
                      else

!         the e-b circular plate element (model=1 and ntype>0)

                         anu21=bx0*ax0/ax1
                         di=(bx1**3)/12.0
                         d11=di*ax0/(1.0-bx0*anu21)
                         d22=d11*(ax1/ax0)
                         d12=bx0*d22
                         do 80 j=1,nn
                         if(item.le.2)then
                    elf(j) = elf(j) + const*sf(j)*fx*x
                 endif
                 do 80 i=1,nn
                 bij = const*gddsf(i)*gddsf(j)*d11*x
                 cij = const*sf(i)*sf(j)*cx*x
                 dij = const*(gddsf(i)*gdsf(j)+gdsf(i)*gddsf(j))*d12
                 eij = const*gdsf(i)*gdsf(j)*d22/x
     80          elk(i,j)=elk(i,j) + bij + cij + dij + eij
              endif
          endif
    100   continue
        else


!         coefficients for the timoshenko beam and circular plate (model=2)
!             full integration for bending coefficients

            do 160 ni=1,ngp
            xi=gauspt(ni,ngp)
            call shape1d(h,ielem,npe,xi)
            const=gj*gauswt(ni,ngp)
            x = 0.0
            do 110 j=1,npe
    110     x = x + sf(j)*elx(j)
            if(ntype.eq.0 .or. ntype.eq.2)then

!         the timoshenko beam element (model=2 and ntype=0 or 2)

               bx=bx0+bx1*x
               cx=cx0+cx1*x
               fx=fx0+fx1*x+fx2*x*x
               jj=1
               do 130 j=1,npe
               if(item.le.2)then
                  elf(jj)=elf(jj)+fx*sf(j)*const
               endif
               ii=1
               do 120 i=1,npe
               cij=sf(i)*sf(j)*const
               bij=gdsf(i)*gdsf(j)*const
               elk(ii,jj)    =elk(ii,jj)    +cx*cij
               elk(ii+1,jj+1)=elk(ii+1,jj+1)+bx*bij
               if(item.ne.0)then
                  elm(ii,jj)    =elm(ii,jj)    +ct0*cij
                  elm(ii+1,jj+1)=elm(ii+1,jj+1)+ct1*cij
               endif
    120        ii=ndf*i+1
    130        jj=ndf*j+1
            else

!          timoshenko circular plate element (model=2 and ntype=1 or 3)
!                       ax0=e1, ax1=e2, bx0=anu12, bx1=h

                   anu21=bx0*ax0/ax1
                   cx=cx0+cx1*x
                   fx=fx0+fx1*x
                   di=(bx1**3)/12.0
                   d11=di*ax0/(1.0-bx0*anu21)
                   d22=d11*(ax1/ax0)
                   d12=bx0*d22
                   jj=1
                   do 150 j=1,npe
                   if(item.le.2)then
                      elf(jj)=elf(jj)+fx*sf(j)*const*x
                   endif
                   ii=1
                   do 140 i=1,npe
                   bij = const*gdsf(i)*gdsf(j)*d11*x
                   cij = const*sf(i)*sf(j)*x
                   dij = const*(gdsf(i)*sf(j)+sf(i)*gdsf(j))*d12
                   eij = const*sf(i)*sf(j)*d22/x
                   elk(ii,jj)    =elk(ii,jj)     + cx*cij
                   elk(ii+1,jj+1)=elk(ii+1,jj+1) + bij + dij + eij
                   if(item.ne.0)then
                      elm(ii,jj)    =elm(ii,jj)    +ct0*cij
                      elm(ii+1,jj+1)=elm(ii+1,jj+1)+ct1*cij
                   endif
        140        ii=ndf*i+1
        150        jj=ndf*j+1
                endif
        160     continue


!          reduced integration is used to evaluate the transverse shear terms

                lgp=ngp-1
                do 230 ni=1,lgp
                xi=gauspt(ni,lgp)

                call shape1d(h,ielem,npe,xi)
                const=gj*gauswt(ni,lgp)

                x = 0.0
                do 170 j=1,npe
        170     x = x + sf(j)*elx(j)
                if(ntype.eq.0 .or. ntype.eq.2)then

!              the timoshenko beam element (model=2 and ntype=0 or 2)
!              ax = gak = ax0 + ax1*x (reduced integration)

               ax=ax0+ax1*x
               jj=1
               do 190 j=1,npe
               ii=1
               do 180 i=1,npe
               b11=gdsf(i)*gdsf(j)*const
               b01=sf(i)*gdsf(j)*const
               b10=gdsf(i)*sf(j)*const
               b00=sf(i)*sf(j)*const
               elk(ii,jj)    =elk(ii,jj)    +ax*b11
               elk(ii,jj+1) =elk(ii,jj+1) +ax*b10
               elk(ii+1,jj) =elk(ii+1,jj) +ax*b01
               elk(ii+1,jj+1)=elk(ii+1,jj+1)+ax*b00
    180        ii=i*ndf+1
    190        jj=j*ndf+1
            else

!         timoshenko circular plate element (model=2 and ntype=1 or 3)
!                    bx1=h, fx2=g13*k (reduced integration)

              a33=bx1*fx2
              jj=1
              do 210 j=1,npe
              ii=1
              do 200 i=1,npe
              bij = const*gdsf(i)*gdsf(j)*x
              cij = const*sf(i)*sf(j)*x
              dij = const*gdsf(i)*sf(j)*x
              dji = const*sf(i)*gdsf(j)*x
              elk(ii,jj)    =elk(ii,jj)     + a33*bij
              elk(ii,jj+1) =elk(ii,jj+1)    + a33*dij
              elk(ii+1,jj) =elk(ii+1,jj)    + a33*dji
              elk(ii+1,jj+1)=elk(ii+1,jj+1) + a33*cij
    200       ii=ndf*i+1
    210       jj=ndf*j+1
          endif
    230   continue
          if(item.eq.0 .and. ntype.gt.1)then
              call timforce(elf,elx,fx0,fx1,fx2,h,ntype,ne,f3,mxelm)
          endif
        endif

          if(item.gt.2)return
             if(item.eq.1 .or. item.eq.2)then

!         equivalent coefficient matrices for time-dependent problems

                if(item .eq. 1)then
               ! if(item .eq. 1)then

!           alfa-family of time approximation for parabolic equations

                        do 250 j=1,nn
                        sum=0.0
                        do 240 i=1,nn
                        sum=sum+(elm(i,j)-a2*elk(i,j))*elu(i)
        240             elk(i,j)=elm(i,j)+a1*elk(i,j)
        250             elf(j)=(a1+a2)*elf(j)+sum
                     else

!           newmark-family of time approximation for hyperbolic equations

                   if(time.eq.0.0)then
                      do 260 j=1,nn
                      do 260 i=1,nn
                      elf(j)=elf(j)-elk(i,j)*elu(i)
        260           elk(i,j)=elm(i,j)
                   else
                      do 270 j=1,nn
                      do 270 i=1,nn
                      elf(j)=elf(j)+elm(i,j)*(a3*elu(i)+a4*elv(i)+a5*ela(i))
        270           elk(i,j)=elk(i,j)+a3*elm(i,j)
                   endif
                endif
            endif
            return
            end


               subroutine constrnt(neq,nhbw,ndf,ncon,icon,vcon,glk,glm,glf, &
                                   trm,mxneq)
!           ____________________________________________________________________

!           the subroutine is called in main to implement specified constraint
!            conditions (e.g., inclined supports) on the condensed system of
!            equations. array glm is used here as a temporary storage array.
!           ____________________________________________________________________

                implicit real*8 (a-h,o-z)
                dimension icon(9),vcon(9),glk(mxneq,mxneq),glf(mxneq), &
                           glm(mxneq,mxneq),trm(mxneq,mxneq)
                pi=3.14159265d0

!           include specified constraint conditions

            do 20 ic=1,neq
               do 10 jc=1,neq
               glm(ic,jc)=0.0
         10    trm(ic,jc)=0.0
    20       trm(ic,ic)=1.0d0
         do 30 ic=1,ncon
             beta=vcon(ic)*pi/180.0d0
             idof=ndf*icon(ic)-1
             trm(idof,idof)    = dcos(beta)
             trm(idof,idof+1) = dsin(beta)
             trm(idof+1,idof) =-dsin(beta)
    30       trm(idof+1,idof+1)= dcos(beta)
             l=0
         do 50 i=1,neq
         do 40 j=1,nhbw
    40   glm(i,l+j)=glk(i,j)
    50   l=l+1
         do 60 i=1,neq
         do 60 j=i,neq
    60   glm(j,i)=glm(i,j)
         do 70 i=1,neq
         do 70 j=1,neq
    70   glk(i,j)=glm(i,j)
         do 80 i=1,neq
         do 80 j=1,neq
         glm(i,j)=0.0
         do 80 k=1,neq
    80   glm(i,j)=glm(i,j)+trm(i,k)*glk(k,j)
         do 90 i=1,neq
         do 90 j=1,neq
         glk(i,j)=0.0
         do 90 k=1,neq
    90   glk(i,j)=glk(i,j)+glm(i,k)*trm(j,k)
         do 100 i=1,neq
         do 100 j=1,neq
   100   trm(i,j)=glk(i,j)
         l=0
         do 120 i=1,neq
         do 110 j=1,nhbw
   110   glk(i,j)=trm(i,l+j)
   120   l=l+1
         do 150 i=1,neq
         glm(i,1)=0.0
         do 140 k=1,neq
   140   glm(i,1)=glm(i,1)+trm(i,k)*glf(k)
   150   glf(i)=glm(i,1)
         return
         end


         subroutine echodata(in,it)
         implicit real*8(a-h,o-z)
         dimension aa(20)
         write(it,40)
        10 continue
           read(in,30,end=20) aa
           write(it,60) aa
           go to 10
        20 continue
           rewind(in)
           write(it,50)
           return
        30 format(20a4)
        40 format(5x,'*** echo of the input data starts ***',/)
        50 format(5x,'**** echo of the input data ends ****',/)
        60 format(1x,20a4)
           end


           subroutine egnsolvr(n,a,b,xx,x,negn,nr,mxneq)
!       __________________________________________________________________

!        the subroutine is called in main to solve the eigenvalue problem

!                                     [a]{x} = lambda[b]{x}

!         the program can be used only for positive-definite [b] matrix.
!         the dimensions of v, vt, w, and ih should be equal to mxneq.
!       __________________________________________________________________

           implicit real*8 (a-h,o-z)
           dimension a(mxneq,mxneq),b(mxneq,mxneq),xx(mxneq),x(mxneq,mxneq)
           dimension v(500,500),vt(500,500),w(500,500),ih(500)

!       call subroutine jacobi to diagonalize [b]

           call jacobi (n,b,negn,nr,v,xx,ih,mxneq)

!       make diagonalized [b] symmetric

           do 10 i=1,n
           do 10 j=1,n
        10 b(j,i)=b(i,j)

!       check (to make sure) that [b] is positive-definite

           do 30 i=1,n
           if (b(i,i))20,30,30
        20 write(6,80)
           stop
        30 continue
!       the eigenvectors of [b] are stored in array v(i,j)
!       form the transpose of [v] as [vt]

       do 40 i=1,n
       do 40 j=1,n
    40 vt(i,j)=v(j,i)

!       find the product [f]=[vt][a][v] and store in [a] to save storage

        call matrxmlt (mxneq,n,vt,a,w)
        call matrxmlt (mxneq,n,w,v,a)

!       get [gi] from diagonalized [b], but store it in [b]

       do 50 i=1,n
    50 b(i,i)=1.0/dsqrt(b(i,i))

!       find the product [q]=[gi][f][gi]=[b][a][b] and store in [a]

        call matrxmlt (mxneq,n,b,a,w)
        call matrxmlt (mxneq,n,w,b,a)

!       we now have the form [q]{z}=lamda{z}. diagonalize [q] to obtain
!       the eigenvalues by calling jacobi. the eigenvalues are returned
!       as diag [a].

       call jacobi (n,a,negn,nr,vt,xx,ih,mxneq)
       do 60 j=1,n
    60 xx(j)=a(j,j)

!       the eigenvectors are computed from the relation,
!                       {x}=[v][gi]{z}=[v][b][vt]
!       since {z} is stored in [vt]

        call matrxmlt (mxneq,n,v,b,w)
        call matrxmlt (mxneq,n,w,vt,x)

    80 format(/'*** matrix [glm] is not positive-definite ***')
       return
       end


        subroutine eqnsolvr(nrm,ncm,neqns,nbw,band,rhs,ires)
!       _________________________________________________________________

!       the subroutine is called in main to solve symmetric and banded set
!       of equations using the gauss elimination method:[band]{u} = {rhs}.
!       the coefficient matrix is input as band(neqns,nbw) and the column
!       vector is input as rhs(neqns), where neqns is the actual number
!       of equations and nbw is the half band width. the true dimensions
!       of the matrix [band] in the calling program, are nrm by ncm. when
!       ires is greater than zero, the right hand elimination is skipped.
!       _________________________________________________________________

             implicit real*8(a-h,o-z)
             dimension band(nrm,ncm),rhs(nrm)

             meqns=neqns-1
             if(ires.le.0) then
                do 30 npiv=1,meqns
                npivot=npiv+1
                lstsub=npiv+nbw-1
                if(lstsub.gt.neqns) then
                   lstsub=neqns
                endif

                do 20 nrow=npivot,lstsub
                ncol=nrow-npiv+1
                factor=band(npiv,ncol)/band(npiv,1)
                do 10 ncol=nrow,lstsub
                icol=ncol-nrow+1
                jcol=ncol-npiv+1
        10      band(nrow,icol)=band(nrow,icol)-factor*band(npiv,jcol)
        20      rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
        30      continue

           else
        40    do 60 npiv=1,meqns
              npivot=npiv+1
              lstsub=npiv+nbw-1
              if(lstsub.gt.neqns) then
                 lstsub=neqns
              endif
              do 50 nrow=npivot,lstsub
              ncol=nrow-npiv+1
              factor=band(npiv,ncol)/band(npiv,1)
        50    rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
        60    continue
           endif

!         back substitution

             do 90 ijk=2,neqns
             npiv=neqns-ijk+2
             rhs(npiv)=rhs(npiv)/band(npiv,1)
             lstsub=npiv-nbw+1
             if(lstsub.lt.1) then
                lstsub=1
             endif
             npivot=npiv-1
       do 80 jki=lstsub,npivot
       nrow=npivot-jki+lstsub
       ncol=npiv-nrow+1
       factor=band(nrow,ncol)
    80 rhs(nrow)=rhs(nrow)-factor*rhs(npiv)
    90 continue
       rhs(1)=rhs(1)/band(1,1)
       return
       end


         subroutine jacobi (n,q,jvec,m,v,x,ih,mxneq)
!        __________________________________________________________________

!          called in egnsolvr to diagonalize [q] by successive rotations

!          description of the variables:

!            n   .... order of the real, symmetric matrix [q] (n > 2)
!           [q] .... the matrix to be diagonalized (destroyed)
!           jvec .... 0, when only eigenvalues alone have to be found
!           [v] .... matrix of eigenvectors
!            m   .... number of rotations performed
!        __________________________________________________________________

         implicit real*8 (a-h,o-z)
         dimension q(mxneq,mxneq),v(mxneq,mxneq),x(mxneq),ih(mxneq)
         epsi=1.0d-08
         if(jvec)10,50,10
    10   do 40 i=1,n
         do 40 j=1,n
         if(i-j)30,20,30
    20   v(i,j)=1.0
         go to 40
    30   v(i,j)=0.0
    40   continue
    50   m=0
         mi=n-1
         do 70 i=1,mi
         x(i)=0.0
         mj=i+1
         do 70 j=mj,n
         if(x(i)-dabs(q(i,j)))60,60,70
    60   x(i)=dabs(q(i,j))
         ih(i)=j
    70   continue
    75   do 100 i=1,mi
         if(i-1)90,90,80
    80   if(xmax-x(i))90,100,100
    90   xmax=x(i)
             ip=i
             jp=ih(i)
        100 continue
             if(xmax-epsi)500,500,110
        110 m=m+1
             if(q(ip,ip)-q(jp,jp))120,130,130
        120 tang=-2.0*q(ip,jp)/(dabs(q(ip,ip)-q(jp,jp))+dsqrt((q(ip,ip) &
                  -q(jp,jp))**2+4.0*q(ip,jp)**2))
             go to 140
        130 tang= 2.0*q(ip,jp)/(dabs(q(ip,ip)-q(jp,jp))+dsqrt((q(ip,ip) &
                  -q(jp,jp))**2+4.0*q(ip,jp)**2))
        140 cosn=1.0/dsqrt(1.0+tang**2)
             sine=tang*cosn
             qii=q(ip,ip)
             q(ip,ip)=cosn**2*(qii+tang*(2.*q(ip,jp)+tang*q(jp,jp)))
             q(jp,jp)=cosn**2*(q(jp,jp)-tang*(2.*q(ip,jp)-tang*qii))
             q(ip,jp)=0.0
             if (q(ip,ip)-q(jp,jp)) 150,190,190
        150 temp=q(ip,ip)
             q(ip,ip)=q(jp,jp)
             q(jp,jp)=temp
             if(sine) 160,170,170
        160 temp=cosn
             goto 180
        170 temp=-cosn
        180 cosn=dabs(sine)
             sine=temp
        190 do 260 i=1,mi
             if (i-ip) 210,260,200
        200 if (i-jp) 210,260,210
        210 if (ih(i)-ip) 220,230,220
        220 if (ih(i)-jp) 260,230,260
        230 k=ih(i)
             temp=q(i,k)
             q(i,k)=0.0
             mj=i+1
             x(i)=0.0
             do 250 j=mj,n
             if (x(i)-dabs(q(i,j))) 240,240,250
        240 x(i)=dabs(q(i,j))
             ih(i)=j
        250 continue
             q(i,k)=temp
        260 continue
             x(ip)=0.0
             x(jp)=0.0
             do 430 i=1,n
             if(i-ip) 270,430,320
        270 temp=q(i,ip)
             q(i,ip)=cosn*temp+sine*q(i,jp)
        if (x(i)-dabs(q(i,ip))) 280,290,290
    280 x(i)=dabs(q(i,ip))
        ih(i)=ip
    290 q(i,jp)=-sine*temp+cosn*q(i,jp)
        if (x(i)-dabs(q(i,jp))) 300,430,430
    300 x(i)=dabs(q(i,jp))
        ih(i)=jp
        go to 430
    320 if(i-jp) 330,430,380
    330 temp=q(ip,i)
        q(ip,i)=cosn*temp+sine*q(i,jp)
        if(x(ip)-dabs(q(ip,i)))340,350,350
    340 x(ip)=dabs(q(ip,i))
        ih(ip)=i
    350 q(i,jp)=-sine*temp+cosn*q(i,jp)
        if (x(i)-dabs(q(i,jp))) 300,430,430
    380 temp=q(ip,i)
        q(ip,i)=cosn*temp+sine*q(jp,i)
        if(x(ip)-dabs(q(ip,i)))390,400,400
    390 x(ip)=dabs(q(ip,i))

        ih(ip)=i
    400 q(jp,i)=-sine*temp+cosn*q(jp,i)
        if(x(jp)-dabs(q(jp,i)))410,430,430
    410 x(jp)=dabs(q(jp,i))
        ih(jp)=i
    430 continue
        if(jvec)440,75,440
    440 do 450 i=1,n
        temp=v(i,ip)
        v(i,ip)=cosn*temp+sine*v(i,jp)
    450 v(i,jp)=-sine*temp+cosn*v(i,jp)
        goto 75
    500 return
        end



        subroutine matrxmlt(mxneq,n,a,b,c)
!       ________________________________________________________________

!       called in egnsolvr to compute the product of matrices [a] & [b]:
!                                   [c]=[a][b]
!       ________________________________________________________________

        implicit real*8 (a-h,o-z)
        dimension a(mxneq,mxneq),b(mxneq,mxneq),c(mxneq,mxneq)
        do 10 i=1,n
        do 10 j=1,n
        c(i,j)=0.0
           do 10 k=1,n
        10 c(i,j)=c(i,j)+a(i,k)*b(k,j)
           return
           end


           subroutine mesh1d(nem,npe,nod,mxelm,mxnod,dx,glx)
!       __________________________________________________________________

!       the subroutine is called in main to compute arrays {glx} and [nod]

!          {glx}.... vector of global coordinates
!          {dx}..... vector of element lengths [dx(1) = node 1 coordinate]
!          [nod].... connectivity matrix
!       __________________________________________________________________

           implicit real*8 (a-h,o-z)
           dimension glx(mxnod),dx(mxnod),nod(mxelm,4)

!       generate the elements of the connectivity matrix

           do 10 i=1,npe
        10 nod(1,i)=i
           do 20 n=2,nem
           do 20 i=1,npe
        20 nod(n,i) = nod(n-1,i)+npe-1

!       generate global coordinates of the global nodes

           glx(1)=dx(1)
           if(npe.eq.2)then
               do 30 i=1,nem
        30     glx(i+1) = glx(i) + dx(i+1)
           else
               do 40 i=1,nem
               ii=2*i
               glx(ii) = glx(ii-1) + 0.5*dx(i+1)
        40     glx(ii+1)=glx(ii-1) + dx(i+1)
           endif
           return
           end

           subroutine postproc(dcax,dcbx,dccx,f3,glf,glx,nod,icont,ielem,npe, &
                               model,ntype,item,mxelm,mxneq,mxnod,nem,ndf)
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

          implicit real*8 (a-h,o-z)
          dimension dcax(mxelm,2),dcbx(mxelm,2),dccx(mxelm,2)
          dimension f3(mxelm),glf(mxneq),glx(mxnod),nod(mxelm,4)
          dimension xp(9),elx(4),elu(9)
          common/io/in,it
          common/shp/sf(4),gdsf(4),gddsf(4),gj
          common/stf2/a1,a2,a3,a4,a5,ax0,ax1,bx0,bx1,cx0,cx1,ct0,ct1,fx0, &
                      fx1,fx2
          data xp/-1.0d0, -0.750d0, -0.50d0, -0.250d0, 0.0d0, 0.250d0, &
                                     0.50d0, 0.750d0, 1.0d0/

         npts=9
         do 80 ne = 1, nem
         if(icont.ne.1) then
             ax0=dcax(ne,1)
             ax1=dcax(ne,2)
             bx0=dcbx(ne,1)
             bx1=dcbx(ne,2)
             cx0=dccx(ne,1)
             cx1=dccx(ne,2)
         endif
         l=0
         do 10 i=1,npe
         ni=nod(ne,i)
         if(icont.ne.1)then
             elx(1)=0.0
             elx(2)=0.5*glx(ne)
             elx(npe)=glx(ne)
         else
             elx(i)=glx(ni)
         endif
         li=(ni-1)*ndf
         do 10 j=1,ndf
           li=li+1
           l=l+1
        10 elu(l)=glf(li)
           h = elx(npe) - elx(1)

           do 70 ni=1,npts
           xi = xp(ni)
           call shape1d(h,ielem,npe,xi)


           if(model.eq.3)then
              w=0.0
              dw=0.0
              ddw=0.0
              xc=elx(1)+0.5*h*(1.0+xi)
              do 20 i=1,4
              w =w + sf(i)*elu(i)
              dw =dw + gdsf(i)*elu(i)
        20    ddw=ddw+ gddsf(i)*elu(i)
              dddw=((elu(1)-elu(3))*2.0/h-(elu(4)+elu(2)))*6.0/(h*h)
              theta=-dw
              if(ntype.eq.0)then
                 bm=-(bx0+xc*bx1)*ddw
                 vf=-(bx0+xc*bx1)*dddw - bx1*ddw
                 write(it,90)xc,w,theta,bm,vf
              else
                 anu21=bx0*ax0/ax1
                 di=(bx1**3)/12.0
                 d11=di*ax0/(1.0-bx0*anu21)
                 d22=d11*(ax1/ax0)
                 d12=bx0*d22
                 bmr=-(d11*ddw*xc+d12*dw)
                 bmt=-(d12*ddw*xc+d22*dw)
                 if(xc.ne.0.0)then
                     sfv=-d11*(xc*dddw+ddw)+d22*dw/xc
                     write(it,90)xc,w,theta,bmr,bmt,sfv
                 else
                     write(it,90)xc,w,theta,bmr,bmt
                 endif
              endif
           else
              xc=0.0
              do 30 i=1,npe
        30    xc=xc+sf(i)*elx(i)
              if(model.eq.1)then
                 u=0.0
                 du=0.0
                 do 40 i=1,npe
                 u=u+sf(i)*elu(i)
        40       du=du+gdsf(i)*elu(i)
               if(ntype.eq.0)then
                  sv=(ax0+ax1*xc)*du
                  write(it,90)xc,u,sv
               else
                  anu21=bx0*ax0/ax1
                  if(ntype.eq.1)then
                      c11=bx1*ax0/(1.0-bx0*anu21)
                      c22=c11*(ax1/ax0)
                      c12=bx0*c22
                  else
                      denom=1.0-bx0-anu21
                      c11=bx1*ax0*(1.0-bx0)/(1.0+bx0)/denom
                      c22=bx1*ax1*(1.0-anu21)/(1.0+anu21)/denom
                      c12=bx0*c22
                  endif
                  if(xc.ne.0.0)then
                      sr=c11*du+c12*u/xc
                      st=c12*du+c22*u/xc
                      write(it,90)xc,u,sr,st
                  else
                     write(it,90)xc,u,du
                  endif
               endif
            else

!       model.eq.2   calculations

                 if(item.eq.0 .and. ntype.gt.1)then
                    h=elx(npe)-elx(1)
                    call timstres(ax0,elu,xi,w,dw,si,dsi,ne,f3,h,mxelm)
                 else
                    w   =0.0
                    dw =0.0
                    psi =0.0
                    dpsi =0.0
                    do 50 i=1,npe
                    l=2*i-1
                    w   = w   +sf(i)*elu(l)
                    dw = dw +gdsf(i)*elu(l)
                    psi = psi +sf(i)*elu(l+1)
    50              dpsi= dpsi+gdsf(i)*elu(l+1)
                 endif
                 if(ntype.eq.0 .or. ntype.eq.2)then
                    bm=(bx0+bx1*xc)*dpsi
                    vf=(ax0+ax1*xc)*(dw+psi)
                    write(it,90)xc,w,psi,bm,vf
                 else
                    anu21=bx0*ax0/ax1
                    di =(bx1**3)/12.0
                    d11=di*ax0/(1.0-bx0*anu21)
                      d22=d11*(ax1/ax0)
                      d12=bx0*d22
                      bmr=(d11*dpsi*xc+d12*psi)
                      bmt=(d12*dpsi*xc+d22*psi)
                      sfv=fx2*(dw+psi)*xc
                      write(it,90)xc,w,psi,bmr,bmt,sfv
                   endif
                endif
            endif
         70 continue
         80 continue
            return
         90 format(2x,6e13.5)
            end

            subroutine reaction(mxelm,mxneq,ndf,nem,nod,npe,ntype,pr,glf, &
                                se,sl,sa,si,cs,sn,cnt,snt,hf,vf,pf,xb)
!        __________________________________________________________________

!         the subroutine is called in main to compute generalized reaction
!        forces in each element of truss (ndf=2) or frame (ndf=3) structure
!        __________________________________________________________________

            implicit real*8(a-h,o-z)
            dimension pr(mxelm),se(mxelm),sl(mxelm),sa(mxelm),si(mxelm)
            dimension cs(mxelm),sn(mxelm),cnt(mxelm),snt(mxelm)
            dimension hf(mxelm),vf(mxelm),pf(mxelm),xb(mxelm)
            dimension nod(mxelm,4),glf(mxneq),elr(6)
            common/stf1/elk(9,9),elm(9,9),elf(9),elx(4),elu(9),elv(9),ela(9)

            nn=npe*ndf
            do 140 n=1,nem
            cn1=cs(n)
            sn1=sn(n)

!        call transfrm to compute element stiffness matrix and force vector

             l=0
             do 100 i=1,npe
             ni=nod(n,i)
             li=(ni-1)*ndf
             do 100 j=1,ndf
             li=li+1
             l=l+1
        100 elu(l)=glf(li)
             call transfrm(mxelm,n,ntype,pr,se,sl,sa,si,cs,sn, &
                           cnt,snt,hf,vf,pf,xb)


!        compute the force and moment resultants

          do 120 i=1,nn
          elr(i) = 0.0
          do 110 j=1,nn
    110   elr(i) = elr(i) + elk(i,j)*elu(j)
    120   elr(i) = elr(i) - elf(i)
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
          write(6,150)n, (elf(i),i=1,nn)
          write(6,160)    (elr(i),i=1,nn)
    140   continue
          return
    150   format (3x,i2,6e12.4)
    160   format (5x,6e12.4,/)
          end



          subroutine shape1d(h,ielem,npe,xi)
!         __________________________________________________________________

!           called in main to compute shape functions and their derivatives
!         for hermite cubic and lagrange linear, quadratic and cubic elements

!            x......... global (i.e., problem) coordinate
!            xi ....... local (i.e., element) coordinate
!            h......... element length
!            {sf}...... interpolation (or shape) functions
!            {dsf}..... first derivative of sf w.r.t. xi
!            {ddsf}.... second derivative of sfh w.r.t. xi
!            {gdsf}.... first derivative of sf w.r.t. x
!            {gddsf}... second derivative of sfh w.r.t. x
!            gj........ determinant of the jacobian matrix
!         __________________________________________________________________

          implicit real*8 (a-h,o-z)
          common/shp/sf(4),gdsf(4),gddsf(4),gj
          dimension dsf(4),ddsf(4)
          if(ielem.eq.0)then

!         hermite interpolation functions (for the euler-bernoulli theory)

            net=4
            sf(1) = 0.25*(2.0-3.0*xi+xi**3)
            sf(2) = -h*(1.0-xi)*(1.0-xi*xi)/8.0
            sf(3) = 0.25*(2.0+3.0*xi-xi**3)
            sf(4) = h*(1.0+xi)*(1.0-xi*xi)/8.0
            dsf(1) = -0.75*(1.0-xi*xi)
            dsf(2) = h*(1.0+2.0*xi-3.0*xi*xi)/8.0
            dsf(3) = 0.75*(1.0-xi*xi)
            dsf(4) = h*(1.0-2.0*xi-3.0*xi*xi)/8.0
            ddsf(1)= 1.5*xi
            ddsf(2)= 0.25*h*(1.0-3.0*xi)
            ddsf(3)= -1.5*xi
            ddsf(4)= -0.25*(1.0+3.0*xi)*h
         else
            net=npe
            if(ielem.eq.1)then


!     lagrange interpolation functions used for linear, quadratic and
!              cubic approximation of second-order equations

!     linear interpolation functions

                sf(1) = 0.5*(1.0-xi)
                sf(2) = 0.5*(1.0+xi)
                dsf(1) = -0.5
                dsf(2) = 0.5
                ddsf(1)= 0.0
                ddsf(2)= 0.0
             else
                if(ielem.eq.2)then

!     quadratic interpolation functions

                     sf(1) =      -0.5*xi*(1.0-xi)
                     sf(2) =      1.0-xi*xi
                     sf(3) =      0.5*xi*(1.0+xi)
                     dsf(1) =     -0.5*(1.0-2.0*xi)
                     dsf(2) =     -2.0*xi
                     dsf(3) =     0.5*(1.0+2.0*xi)
                     ddsf(1)=     1.0
                     ddsf(2)=     -2.0
                     ddsf(3)=     1.0
                  else

!     cubic interpolation functions

                      sf(1)    = 0.0625*(1.0-xi)*(9.0*xi*xi-1.)
                      sf(2)    = 0.5625*(1.0-xi*xi)*(1.0-3.0*xi)
                  sf(3) =    0.5625*(1.0-xi*xi)*(1.0+3.0*xi)
                  sf(4) =    0.0625*(9.0*xi*xi-1.0)*(1.0+xi)
                  dsf(1) =   0.0625*(1.0+18.0*xi-27.0*xi*xi)
                  dsf(2) =   0.5625*(-3.0-2.0*xi+9.0*xi*xi)
                  dsf(3) =   0.5625*(3.0-2.0*xi-9.0*xi*xi)
                  dsf(4) =   0.0625*(18.0*xi+27.0*xi*xi-1.0)
                  ddsf(1)=   0.0625*(18.0-54.0*xi)
                  ddsf(2)=   0.5625*(-2.0+18.0*xi)
                  ddsf(3)=   0.5625*(-2.0-18.0*xi)
                  ddsf(4)=   0.0625*(18.0+54.0*xi)
               endif
            endif
         endif

!        compute derivatives of the interpolation functions w.r.t. x

     80 gj = h*0.5
        do 90 i = 1,net
        gdsf(i) = dsf(i)/gj
     90 gddsf(i) = ddsf(i)/gj/gj
        return
        end


         subroutine timforce(elf,elx,fx0,fx1,fx2,h,ntype,ne,f3,mxelm)
!        __________________________________________________________________

!           called in coeffcnt to compute element force vector for the
!                consistent interpolation timoshenko element (cie)
!        __________________________________________________________________

         implicit real*8(a-h,o-z)
         common/shp/sf(4),gdsf(4),gddsf(4),gj
         dimension gauspt(5,5),gauswt(5,5),elf(9),elx(4),ex(3),f3(mxelm)

         data gauspt/5*0.0d0,-.57735027d0,.57735027d0,3*0.0d0,-.77459667d0, &
          0.0d0,.77459667d0,2*0.0d0,-.86113631d0,-.33998104d0,.33998104d0, &
         .86113631d0,0.0d0,-.906180d0,-.538469d0,0.0d0,.538469d0,.906180d0/

         data gauswt/2.0d0,4*0.0d0,2*1.0d0,3*0.0d0,.55555555d0,.88888888d0, &
          0.55555555d0,2*0.0d0,.34785485d0,2*.65214515d0,.34785485d0,0.0d0, &
          0.236927d0,.478629d0,.568889d0,.478629d0,.236927d0/


        npe=3
        iel=2
        ndf=2
        ngp=iel+1
        do 10 i=1,6
10      elf(i)=0.0

         ex(1)=elx(1)
         ex(2)=elx(1)+0.5*h
         ex(3)=elx(2)

         do 50 ni=1,ngp
         xi=gauspt(ni,ngp)
         call shape1d(h,iel,npe,xi)
         const=gj*gauswt(ni,ngp)
         x = 0.0
         do 20 j=1,npe
   20    x = x + sf(j)*ex(j)

!     compute the polynomial variation of fx

         if(ntype.eq.2)then
            fx=fx0+(fx1+fx2*x)*x
         else
            fx=(fx0+fx1*x)*x
         endif

!     element force vector for the consistent interpolation beam element

   25    ii=1
         do 40 i=1,npe
         elf(ii)=elf(ii)+fx*sf(i)*const
   40    ii=ndf*i+1
   50    continue

!     rearrange the element coefficients

         f3(ne)=elf(3)
         elf(1)=elf(1)+0.5*f3(ne)
         elf(2)=-0.125*f3(ne)*h
         elf(3)=elf(5)+0.5*f3(ne)
         elf(4)= 0.125*f3(ne)*h
         return
         end


         subroutine timstres(ga,elu,xi,w,dw,s,ds,ne,f3,h,mxelm)
!     __________________________________________________________________

!     called in postproc to compute solution and its global derivatives
!       at nine points (including the nodes) of the timoshenko element

!         xc........     global (i.e., problem) coordinate
!         xi .......     local (i.e., element) coordinate
!         sfl, sfq..     lagrange linear and quadratic shape functions
!         dsfl,dsfq:     first derivative of sf w.r.t. global coordinate
!              elu....... column vector of generalized displacements
!              w, dw..... transverse deflection and its derivative
!              s, ds..... rotation and its derivative
!           __________________________________________________________________

         implicit real*8 (a-h,o-z)
         common/io/in,it
         dimension elu(9),sfl(2),sfq(3),dsfl(2),dsfq(3),f3(mxelm)

         gj =       h*0.5



!           interpolation functions for the lagrange linear element

         sfl(1) = 0.5*(1.0-xi)
         sfl(2) = 0.5*(1.0+xi)
         dsfl(1) = -0.5/gj
         dsfl(2) = 0.5/gj

!           interpolation functions for the lagrange quadratic element

         sfq(1) = -0.5*xi*(1.0-xi)
         sfq(2) = 1.0-xi*xi
         sfq(3) = 0.5*xi*(1.0+xi)
         dsfq(1) = -0.5*(1.0-2.0*xi)/gj
         dsfq(2) = -2.0*xi/gj
         dsfq(3) = 0.5*(1.0+2.0*xi)/gj

         w3=(3.0*h*f3(ne)/ga + 8.0*(elu(1)+elu(3)) &
                             + 2.0*(elu(4)-elu(2))*h)/16.0
         w = sfq(1)*elu(1) + sfq(2)*w3 + sfq(3)*elu(3)
         dw= dsfq(1)*elu(1) +dsfq(2)*w3 +dsfq(3)*elu(3)
         s = sfl(1)*elu(2) + sfl(2)*elu(4)
         ds= dsfl(1)*elu(2) +dsfl(2)*elu(4)

         return
         end


         subroutine transfrm(mxelm,n,ntype,pr,se,sl,sa,si,cs,sn,cnt,snt, &
                             hf,vf,pf,xb)
!        __________________________________________________________________

!           called in both main and reaction to compute stiffness matrix and
!            force vector for the truss (ndf=2) and frame (ndf=3) elements

!           se......young's modulus
!           sl......element length
!           sa......cross-sectional area
!         si......moment of inertia
!         cs......cosine of the angle of orientation
!         sn......sine of the angle of orientation
!         hf......distributed force along the length of the element
!         vf......distributed force transverse to the element
!         pf......point force at point other than nodes
!         xb......distance along the length from node 1 of the element
!                 of the location of the point force, pf
!         cnt,snt:direction cosines of the point force's line of application
!         __________________________________________________________________

             implicit real*8(a-h,o-z)
             dimension pr(mxelm),se(mxelm),sl(mxelm),sa(mxelm),si(mxelm)
             dimension cs(mxelm),sn(mxelm),cnt(mxelm),snt(mxelm)
             dimension hf(mxelm),vf(mxelm),pf(mxelm),xb(mxelm)
             dimension trm(6,6),tmpk(6,6)
             common/stf1/elk(9,9),elm(9,9),elf(9),elx(4),elu(9),elv(9),ela(9)

             cn1=cs(n)
             sn1=sn(n)
             cn2=cn1*cn1
             sn2=sn1*sn1
             csn=cn1*sn1

!         element coefficients

             if(ntype.eq.0) then



!         the plane truss element

                nn=4
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

                do 10 i=1,nn
                do 10 j=i,nn
        10      elk(i,j) = elk(j,i)

!         contribution of the point force to nodal forces

            xi=xb(n)/sl(n)
            sfl1 = 1.0-xi
            sfl2 = xi

            f1=0.5*hf(n)*sl(n)
            f3=0.5*hf(n)*sl(n)
            elf(1) = f1*cn1
            elf(2) = f1*sn1
            elf(3) = f3*cn1
            elf(4) = f3*sn1
         else
            nn=6
            if(ntype.eq.1)then

!        the euler-bernoulli frame element

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

                 do 20 i=1,nn
                 do 20 j=i,nn
    20           elk(i,j) = elk(j,i)

!     contribution of the point force to nodal generalized forces

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
             else

!     the timoshenko frame element (shear coefficient=5/6)

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

                 do 25 i=1,nn
                 do 25 j=1,nn
    25           trm(j,i)=0.0

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


                 do 30 i=1,nn
                 do 30 j=i,nn
    30           elk(i,j) = elk(j,i)

                 do 40 i=1,nn
                 do 40 j=1,nn
                 tmpk(i,j)=0.0
                 do 40 k=1,nn
    40           tmpk(i,j)=tmpk(i,j)+trm(k,i)*elk(k,j)

                 do 50 i=1,nn
                 do 50 j=1,nn
                 elk(i,j)=0.0
                 do 50 k=1,nn
    50           elk(i,j)=elk(i,j)+tmpk(i,k)*trm(k,j)

!        contribution of the point force to nodal generalized forces

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
             endif
         endif
         return
end
