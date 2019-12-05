! Uses Lambeck et al. 2014 sea-level record as forcing for LGM simulations
! Written by Clemens Schannwell
! Last Updated 26.11.2018 added dynamic memory allocation (awk format is
! %9.3f%9.3f)
FUNCTION GetSeaLevel ( Model, nodenumber, Time) RESULT(SeaLevelOut)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   INTEGER :: nodenumber, j,k, indIt, io, i, LineCount
   REAL(KIND=dp) :: Time, SeaLevelOut
   !REAL, DIMENSION(326) :: SLTimeIn, SLCIn, SLTime
   REAL, DIMENSION(:), ALLOCATABLE :: SLTimeIn, SLCIn, SLTime
   LOGICAL :: FirstTime = .TRUE.

   SAVE :: FirstTime, SLTime, SLCIn
!------------------------------------------------------------------------------
! Real function starts
!------------------------------------------------------------------------------
! Read data from file
!------------------------------------------------------------------------------
  IF (FirstTime) THEN
          FirstTime = .FALSE.
          ! Dry read to get dimension for dynamic memory allocation
          OPEN(11,FILE="src/SeaLevel/SeaLevelGlacialCycle",STATUS='old',ACTION='read')
          LineCount = 0
          DO
            READ(11,*,iostat=io)
            IF (io/=0) EXIT
            LineCount = LineCount + 1
          END DO
          REWIND(11)
          ALLOCATE (SLTimeIn(LineCount), SLCIn(LineCount), SLTime(LineCount))
          DO i=1,LineCount
            READ(11,'(F9.3,F9.3)') SLTimeIn(i), SLCIn(i)
            !print *, 'TIMEVEC ', SLTimeIn(i), SLCIn(i)
          END DO
          DO j = 1, size(SLTimeIn)
              SLTime(j) = SLTimeIn(j)*1000
          END DO
  END IF
  indIt=0
  ! Linear interpolation onto exact simulation time
  DO WHILE (Time > SLTime(indIt))
      indIt = indIt + 1
  END DO
  SeaLevelOut = (Time - SLTime(indIt-1)) / (SLTime(indIt) - SLTime(indIt-1)) &
          * (SLCIn(indIt)-SLCIn(indIt-1)) + SLCIn(indIt-1)
  !PRINT *, 'PRINT: ', SeaLevelOut

END FUNCTION GetSeaLevel
