! Get BMB as Function of Ice Shelf thickness and ocean temperature
! This follows Goose & Beckmann 2003; implementation follows Martin et al. 2011
!(TC) or Dolan 2018 (NatComm)
! Written by Clemens Schannwell
! Last Updated 26.11.2018 added dynamic memory allocation (awk format is
! %12.5f%12.5f)
FUNCTION GetBMB ( Model, nodenumber, Time) RESULT(BMBOut)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: DepthSol
   TYPE(ValueList_t), POINTER :: Material
   INTEGER, POINTER :: DepthPerm(:)
   REAL(kind=dp), POINTER :: DepthVal(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber, indIt, io, i, LineCount
   REAL(KIND=dp) :: Time, BMBOut, SO, Tf, rhoi, rhow, cpo, yt, L, Fmelt
   REAL(KIND=dp) :: TocInit, DeltaTOcean, Toc
   !REAL, DIMENSION(2941) :: OceanTempTimeIn, OceanTempIn
   REAL, DIMENSION(:), ALLOCATABLE :: OceanTempTimeIn, OceanTempIn
   REAL(KIND=dp), ALLOCATABLE :: BMB0(:)
   LOGICAL :: FirstTime=.True., UnFoundFatal, GotIt

   SAVE :: FirstTime, OceanTempIn, OceanTempTimeIn
!------------------------------------------------------------------------------
! Real function starts
!------------------------------------------------------------------------------
! Read data from file
!------------------------------------------------------------------------------
  !Time = Time +15010
  IF (FirstTime) THEN
          FirstTime = .FALSE.
          ! Dry read to get dimension for dynamic memory allocation
          OPEN(11,FILE="src/BMB/BMBForcingTimes2GlacialCycle",STATUS='old',ACTION='read')
          LineCount = 0
          DO
            READ(11,*,iostat=io)
            IF (io/=0) EXIT
            LineCount = LineCount + 1
          END DO
          REWIND(11)
          ALLOCATE (OceanTempTimeIn(LineCount), OceanTempIn(LineCount))
          DO i=1,LineCount
            READ(11,'(F12.5,F12.5)') OceanTempTimeIn(i), OceanTempIn(i)
          END DO
  END IF
!------------------------------------------------------------------------------
! Get all constants for calculation of BMB
!------------------------------------------------------------------------------
  rhow = GetConstReal( Model % Constants, 'Water Density', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No water density is specified in sif file.')
  END IF
  SO = GetConstReal( Model % Constants, 'Ocean Salinity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No ocean salinity is specified in the constant section of the sif file.')
  END IF
  cpo = GetConstReal( Model % Constants, 'Heat Capacity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No heat capacity is specified in the constant section of the sif file.')
  END IF
  yt = GetConstReal( Model % Constants, 'Thermal Exchange Velocity', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No thermal exchange velocity is specified in the constant section of the sif file.')
  END IF
  L = GetConstReal( Model % Constants, 'Latent Heat of Fusion', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No latent heat of fusion is specified in the constant section of the sif file.')
  END IF
  Fmelt = GetConstReal( Model % Constants, 'FMelt', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No tuning factor is specified in the constant section of the sif file.')
  END IF
  TocInit = GetConstReal( Model % Constants, 'OceanTempInit', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No initial ocean temperature is specified in the constant section of the sif file.')
  END IF
  Material => GetMaterial()
  rhoi = GetConstReal( Material, 'Density', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF BMBBeckmannGoose', 'No ice density is specified in the material section of the sif file.')
  END IF
!------------------------------------------------------------------------------
! Done reading constants
!------------------------------------------------------------------------------
  indIt=0
  DO WHILE (Time > OceanTempTimeIn(indIt))
      indIt = indIt + 1
  END DO
  ! Linear interpolation onto exact simulation time
  DeltaTOcean = (Time - OceanTempTimeIn(indIt-1)) / (OceanTempTimeIn(indIt) &
          - OceanTempTimeIn(indIt-1))* (OceanTempIn(indIt)-OceanTempIn(indIt-1)) &
          + OceanTempIn(indIt-1)
  Toc = TocInit + DeltaTOcean
!------------------------------------------------------------------------------
!Get Depthwithout catching any error messages if fields don't exist
!------------------------------------------------------------------------------
   DepthSol => VariableGet( Model % Variables, 'Depth',UnFoundFatal=UnFoundFatal)
   DepthPerm => DepthSol % Perm
   DepthVal => DepthSol % Values

   ! Compute local freezing temperature
   Tf = 0.0939-0.057*SO + 7.64E-4*-DepthVal(DepthPerm(nodenumber))
   ! Compute basal mass balance
   !BMBOut = (rhow*cpo*yt*Fmelt*((Toc-Tf))/(L*rhoi))*365.25*24*60*60
   BMBOut = (rhow*cpo*yt*Fmelt*(ABS(Toc-Tf)*(Toc-Tf))/(L*rhoi))*365.25*24*60*60
   !print *, 'DELTAO ', DeltaTOcean
   !write(*,*) 'BMBOut:',Model % Nodes % x(nodenumber),Model % Nodes % y(nodenumber),BMBOut

END FUNCTION GetBMB




