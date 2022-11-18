PROGRAM plane_strain_elastic_FE
!-------------------------------------------------------------------------
! Program General two- dimensional analysis of elastic solids in plane
! strain.
!
!-------------------------------------------------------------------------
USE CW3module  ! use library of subroutines and functions
IMPLICIT NONE
!!**********************Initialisation**********************************!!
!-----------------------fixed variables---------------------------------
INTEGER::              &
         i,k,          & ! i,k:simple counter for loops, etc
         iel,          & ! counter for elements
         loaded_nodes, & ! number of loaded nodes
         ndim,         & ! number of dimensions
         ndof,         & ! number of freedoms per element
         nels,         & ! number of elements
         neq,          & ! total number of (non-zero) freedoms in problem
         nip,          & ! number of gauss points
         nn,           & ! number of total nodes in problem
         nod,          & ! number of nodes per element
         nodof,        & ! number of freedoms per node
         nprops=2,     & ! number of material properties
         np_types,     & ! number of different properties
         nr,           & ! number of restrained nodes
         nst,          & ! number of stress(strain) terms
         nlen            ! number of characters in character string
REAL*8::               &
         det,          & ! determinant of the Jacobian matrix
         zero=0.d0
CHARACTER(len=15)::    &
         element,      & ! types of element
         argv            !file name
!-----------------------dynamic arrays------------------------------------
INTEGER,ALLOCATABLE::  &
         etype(:),     & ! element property type vector - i.e. vector with which material type each element is.
         g(:),         & ! element steering vector - i.e. which freedoms belong to which element
         g_g(:,:),     & ! global element steering matrix - large vector of g with all elements freedoms
         g_num(:,:),   & ! global element node numbers matrix
         kdiag(:),     & ! diagonal term location vector
         nf(:,:),      & ! nodal freedom matrix
         node(:),      & ! loaded nodes vector
         num(:)          ! element node numbers vector
REAL*8,ALLOCATABLE::   &
         bee(:,:),     & ! strain-displacement matrix
         coord(:,:),   & ! element nodal coordinates
         dee(:,:),     & ! stress-strain matrix
         der(:,:),     & ! shape function derivatives with respect to local coordinates
         deriv(:,:),   & ! shape function derivatives with respect to global coordinates
         eld(:),       & ! element nodal displacements
         fun(:),       & ! shape functions
         gc(:),        & ! integrating point coordinates
         g_coord(:,:), & ! global nodal coordinates
         jac(:,:),     & ! Jacobian matrix
         km(:,:),      & ! element stiffness matrix
         kv(:),        & ! global stiffness matrix
         loads(:),     & ! nodal loads and displacements
         points(:,:),  & ! integrating point local coordinates
         prop(:,:),    & ! element properties (E and v for each element)
         sigma(:),     & ! stress terms
         weights(:),   & ! weighting coefficients
         val_load(:,:)   ! applied nodal load weightings
!-----------------------input and initialisation--------------------------

CALL getname(argv,nlen)
OPEN(10,FILE=argv(1:nlen)//'.dat')
OPEN(11,FILE=argv(1:nlen)//'.res')

READ(10,*)element,nod,nels,nn,nip,nodof,nst,ndim,np_types
ndof=nod*nodof
ALLOCATE(                                                                 &
   nf(nodof,nn),points(nip,ndim),dee(nst,nst),g_coord(ndim,nn),           &
   coord(nod,ndim),jac(ndim,ndim),weights(nip),num(nod),g_num(nod,nels),  &
   der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),km(ndof,ndof),eld(ndof),   &
   sigma(nst),g(ndof),g_g(ndof,nels),gc(ndim),fun(nod),etype(nels),       &
   prop(nprops,np_types)  )
READ(10,*)prop
etype=1
IF (np_types>1) READ(10,*)etype
READ(10,*)g_coord
READ(10,*)g_num			                  ! global numbering
!
nf=1
READ(10,*)nr,(k,nf(:,k),i=1,nr)
CALL formnf(nf)                           ! to set the boundary conditions
neq=MAXVAL(nf)
ALLOCATE(kdiag(neq),loads(0:neq))
!
!-----------------loop over the elements to find global arrays sizes-----
kdiag=0
DO iel=1,nels
   num=g_num(:,iel)
   CALL num_to_g(num,nf,g)
   g_g(:,iel)=g
   CALL fkdiag(kdiag,g)
END DO elements_1
DO i=2,neq
   kdiag(i)=kdiag(i)+kdiag(i-1)
END DO
ALLOCATE(kv(kdiag(neq)))
WRITE(11,'(2(A,I5))')                                                    &
   " There are",neq," equations and the skyline storage is",kdiag(neq)
!
!!*******************Caluclate stiffness matrix*************************!!
!----------------------element stiffness integration and assembly-------
!
!   <<<   PROGRAM SECTION MISSING --  WRITE ME   >>>
!
!----------------------------loading condition---------------------------
loads=zero
READ (10,*)loaded_nodes
ALLOCATE (node(loaded_nodes),val_load(loaded_nodes,ndim))
READ (10,*)(node(i),val_load(i,:),i=1,loaded_nodes)
DO i=1,loaded_nodes
   loads(nf(:,node(i)))=val_load(i,:)
END DO
!
!!********************************Solve***********************************!!
!-----------------------equation solution---------------------------------
CALL sparin(kv,kdiag)
CALL spabac(kv,loads,kdiag)
!
!!*****************************Output results******************************!!
!---------------------------------Output displacements--------------
!
!   <<<   PROGRAM SECTION MISSING --  WRITE ME   >>>
!
!-----------recover stresses at element Gauss-points and output--------------
!
!   <<<   PROGRAM SECTION MISSING --  WRITE ME   >>>
!

!--------------------output images-----------------------------------------
CALL mesh(g_coord,g_num,argv,nlen,12)
CALL dismsh(loads,nf,0.05d0,g_coord,g_num,argv,nlen,13)
CALL vecmsh(loads,nf,0.05d0,0.1d0,g_coord,g_num,argv,nlen,14)

END PROGRAM plane_strain_elastic_FE

