!      ----------------------------------------------------------------
!                        f   o   r   m   a    t  s
!      ----------------------------------------------------------------
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


350 format(/,8x, 'element type (0, hermite,>0, lagrange)..=',i4,/, &
              8x, 'no. of deg. of freedom per node, ndf....=',i4,/, &
              8x, 'no. of elements in the mesh, nem........=',i4,/, &
              8x, 'no. of total dof in the model, neq......=',i4,/, &
              8x, 'half bandwidth of matrix [glk], nhbw ...=',i4,/, &
              8x, 'no. of specified primary dof, nspv......=',i4,/, &
              8x, 'no. of specified secondary dof, nssv....=',i4,/, &
              8x, 'no. of specified newton b. c.: nnbc.....=',i4,/, &
              8x, 'no. of speci. multi-pt. cond.: nmpc.....=',i4)


390 format(/,8x,'coefficient, ct0........................=',e12.4,/, &
             8x,'coefficient, ct1........................=',e12.4,/, &
             8x,'parameter, alfa.........................=',e12.4,/, &
             8x,'parameter, gama.........................=',e12.4,/, &
             8x,'time increment, dt......................=',e12.4,/, &
             8x,'no. of time steps, ntime................=',i4,/, &
             8x,'time-step interval to print soln., intvl=',i4,/)


400 format(/,8x,'coefficient, ct0........................=',e12.4,/, &
             8x,'coefficient, ct1........................=',e12.4,/)



430 format(/,5x,'properties of element =',i3,//, &
             8x,'element length, h ....... =',e12.4)



440 format( 8x,'ax0 =',e12.4,5x,'ax1 =',e12.4,/, &
             8x,'bx0 =',e12.4,5x,'bx1 =',e12.4,/, &
             8x,'cx0 =',e12.4,5x,'cx1 =',e12.4,/, &
             8x,'fx0 =',e12.4,5x,'fx1 =',e12.4,5x,'fx2 =',e12.4,/)

445 format( 8x,'ax0 =',e12.4,5x,'ax1 =',e12.4,/,  &
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


470 format(8x,'modulus of elasticity,       se........ =',e12.4,/, &
           8x,'length of the element,       sl........ =',e12.4,/, &
           8x,'area of cross section,       sa........ =',e12.4,/, &
           8x,'cosine of orientation,       cn........ =',e12.4,/, &
           8x,'sine of orientation,         sn........ =',e12.4,/, &
           8x,'axial body force (constant), hf........ =',e12.4,/, &
           8x,'nodal connectivity:          nod(i,j).. =',2i6,/)


630 format(/,9x,'generalized forces in the element coordinates',/, &
     5x,'(second line gives the results in the global coordinates)')

640 format(/,3x,'*** the solution has reached a steady state ***', &
            /,3x,'solution at the two consecutive time steps follows:')


660 format(7x,' x ',6x, 'deflect.',5x,'rotation',4x,'moment, mr', &
           3x,'moment, mt',3x,'shear force')
