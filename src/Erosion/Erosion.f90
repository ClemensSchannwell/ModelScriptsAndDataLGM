FUNCTION ComputeErosion ( Model, nodenumber, x) RESULT(Erosion)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: y , x
   INTEGER :: nodenumber

   TYPE(Variable_t), POINTER :: NormalVar, FlowVariable, GroundedMaskVar
   REAL(KIND=dp), POINTER :: NormalValues(:), FlowValues(:), GroundedMaskValues(:)
   INTEGER, POINTER :: NormalPerm(:), FlowPerm(:), GroundedMaskPerm(:)
   INTEGER :: DIM, i
   REAL (KIND=dp) :: Erosion, ErosionExp, ErodabConst
   REAL (KIND=dp) :: ut, un, GLMask
   REAL (KIND=dp), ALLOCATABLE :: normal(:), velo(:)
   LOGICAL :: GotIt, FirstTime = .TRUE., SSA = .FALSE., UnFoundFatal=.TRUE.

   CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

   SAVE :: normal, velo, DIM, SSA
   SAVE :: FlowSolverName, FirstTime

   ErosionExp = GetConstReal( Model % Constants, 'Erosion Exponent', GotIt )
   IF (.NOT.GotIt) THEN
     CALL FATAL('USF Erosion', 'Erosion exponent not found in constant section of sif file.')
   END IF
   ErodabConst = GetConstReal( Model % Constants, 'Erodability Constant', GotIt )
   IF (.NOT.GotIt) THEN
     CALL FATAL('USF Erosion', 'Erodiability constant not found in constant section of sif file.')
   END IF
   !PRINT *, 'ErosionExp', ErosionExp
   IF (FirstTime) THEN
      FirstTime = .FALSE.
      DIM = CoordinateSystemDimension()
      !n = Model % MaxElementNodes
      IF ((DIM == 2).OR.(DIM == 3))  THEN
         ALLOCATE(normal(DIM), velo(DIM))
      ELSE
         CALL FATAL('USF Erosion', 'Bad dimension of the problem')
      END IF

      FlowSolverName = GetString( Model % Solver % Values , 'Flow Solver Name', GotIt )
      IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
      SELECT CASE (FlowSolverName)
      CASE ('ssabasalflow')
         SSA = .TRUE.
      END SELECT
      write(*,*)FlowSolverName, SSA
   END IF
  !n = GetElementNOFNodes()
  ! GroundedMask import
  GroundedMaskVar => VariableGet( Model % Mesh % Variables, 'GroundedMask',UnFoundFatal=UnFoundFatal)
  GroundedMaskValues => GroundedMaskVar % Values
  GroundedMaskPerm => GroundedMaskVar % Perm

  ! Get the variables to compute ut
  FlowVariable => VariableGet( Model % Variables, FlowSolverName,UnFoundFatal=UnFoundFatal)
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values

  ! NS, AIFlow cases
  IF (.NOT.SSA) THEN
     ! Get the variable to compute the normal
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
     NormalPerm => NormalVar % Perm
     NormalValues => NormalVar % Values
   DO i=1, DIM
        normal(i) = -NormalValues(DIM*(NormalPerm(nodenumber)-1) + i)
        velo(i) = FlowValues( (DIM+1)*(FlowPerm(nodenumber)-1) + i )
     END DO
     un = SUM(velo(1:DIM)*normal(1:DIM))
     ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )
     ! SSA Flow case
  ELSE
     DO i=1, DIM-1
        velo(i) = FlowValues( (DIM-1)*(FlowPerm(nodenumber)-1) + i )
     END DO
     ut = SQRT(SUM( velo(1:DIM-1)**2.0 ))

  END IF
  IF (GroundedMaskValues(GroundedMaskPerm(nodenumber)).GT.-1.0_dp) THEN
    Erosion = ErodabConst * ( ut**ErosionExp )
  ELSE
    Erosion = 0.0_dp
  END IF
  !PRINT *, 'GLMask', GroundedMaskValues(GroundedMaskPerm(nodenumber))

END FUNCTION ComputeErosion
