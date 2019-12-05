! Computes the surface temperature assuming a constant temperature increase for
! the region (from ice core data) and a height correction due to changes in the
! geometry
! Written by Clemens Schannwell
! Last Updated 26.11.2018 added dynamic memory allocation (awk format is
! %10.3f%9.3f)
FUNCTION GetSurfTemp ( Model, nodenumber, Time) RESULT(TSurf)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Variable_t), POINTER :: TInitSol,ZsSol,ZsInitSol
   INTEGER :: nodenumber, j,k, indIt, io, i, LineCount
   INTEGER, POINTER :: TInitPerm(:), ZsPerm(:), ZsInitPerm(:)
   REAL(kind=dp), POINTER :: TInitVal(:), ZsVal(:), ZsInitVal(:)
   REAL(KIND=dp) :: Time, DeltaTSurf, DeltaSurf, TSurf, LapseRate
   REAL, DIMENSION(:), ALLOCATABLE :: SurfTempTimeIn, SurfTempIn
   LOGICAL :: FirstTime = .TRUE., UnFoundFatal, GotIt

   SAVE :: FirstTime, SurfTempIn, SurfTempTimeIn
!------------------------------------------------------------------------------
! Real function starts
!------------------------------------------------------------------------------
! Read data from file
!------------------------------------------------------------------------------
  IF (FirstTime) THEN
          FirstTime = .FALSE.
          ! Dry read to get dimension for dynamic memory allocation
          OPEN(11,FILE="src/SurfTemp/SurfTempForcingGlacialCycle",STATUS='old',ACTION='read')
          LineCount = 0
          DO
            READ(11,*,iostat=io)
            IF (io/=0) EXIT
            LineCount = LineCount + 1
          END DO
          REWIND(11)
          ALLOCATE (SurfTempTimeIn(LineCount), SurfTempIn(LineCount))
          DO i=1,LineCount
            READ(11,'(F10.3,F9.3)') SurfTempTimeIn(i), SurfTempIn(i)
            !print *, 'TIMEVEC ', SurfTempTimeIn(i), SurfTempIn(i)
          END DO
  END IF

  indIt=0
  DO WHILE (Time > SurfTempTimeIn(indIt))
      indIt = indIt + 1
  END DO
  DeltaTSurf = (Time - SurfTempTimeIn(indIt-1)) / (SurfTempTimeIn(indIt) &
          - SurfTempTimeIn(indIt-1))* (SurfTempIn(indIt)-SurfTempIn(indIt-1)) &
          + SurfTempIn(indIt-1)
  !Get the initial Temperature profile
  !Get without catching any error messages if fields don't exist
  TInitSol => VariableGet( Model % Variables, 'TempInit',UnFoundFatal=UnFoundFatal)
  TInitPerm => TInitSol % Perm
  TInitVal => TInitSol % Values
  !Get the initial ice surface profile
  !Get without catching any error messages if fields don't exist
  ZsInitSol => VariableGet( Model % Variables, 'ZsInit',UnFoundFatal=UnFoundFatal)
  ZsInitPerm => ZsInitSol % Perm
  ZsInitVal => ZsInitSol % Values
  !Get the current ice surface profile
  !Get without catching any error messages if fields don't exist
  ZsSol => VariableGet( Model % Variables, 'Zs',UnFoundFatal=UnFoundFatal)
  ZsPerm => ZsSol % Perm
  ZsVal => ZsSol % Values
  !Compute difference in ice surface elevation
  DeltaSurf = ZsVal(ZsPerm(nodenumber)) - ZsInitVal(ZsInitPerm(nodenumber))
  !Fetch atmospheric lapse rate from sif file
  LapseRate = GetConstReal( Model % Constants, 'Atmospheric Lapse Rate', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF SurfaceTemp', 'No atmospheric lapse rate is specified in the constant section of the sif file.')
  END IF
  !Compute new surface temperature distribution
  TSurf = TInitVal(TInitPerm(nodenumber)) - LapseRate * DeltaSurf + DeltaTSurf
END FUNCTION GetSurfTemp
