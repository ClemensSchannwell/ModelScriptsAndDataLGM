FUNCTION GetSMB ( Model, nodenumber, x) RESULT(SMB)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Variable_t), POINTER :: TInitSol, TSurfSol, SMBInitSol
   INTEGER :: nodenumber
   INTEGER, POINTER :: TInitPerm(:), TSurfPerm(:), SMBInitPerm(:)
   REAL(kind=dp), POINTER :: TInitVal(:), TSurfVal(:), SMBInitVal(:)
   REAL(KIND=dp) :: DeltaAcc, x, SMB
   LOGICAL :: UnFoundFatal, GotIt

  !Get the initial Temperature profile
  !Get without catching any error messages if fields don't exist
  TInitSol => VariableGet( Model % Variables, 'TempInit',UnFoundFatal=UnFoundFatal)
  TInitPerm => TInitSol % Perm
  TInitVal => TInitSol % Values
  !Get the initial Temperature profile
  !Get without catching any error messages if fields don't exist
  TSurfSol => VariableGet( Model % Variables, 'TSurf',UnFoundFatal=UnFoundFatal)
  TSurfPerm => TSurfSol % Perm
  TSurfVal => TSurfSol % Values
  !Get the initial SMB profile
  !Get without catching any error messages if fields don't exist
  SMBInitSol => VariableGet( Model % Variables, 'SMBInit',UnFoundFatal=UnFoundFatal)
  SMBInitPerm => SMBInitSol % Perm
  SMBInitVal => SMBInitSol % Values
  !!Fetch atmospheric lapse rate from sif file
  DeltaAcc = GetConstReal( Model % Constants, 'Delta Accumulation', GotIt )
  IF (.NOT.GotIt) THEN
     CALL FATAL('USF SMB', 'No delta accumulation is specified in the constant section of the sif file.')
  END IF
  !Compute new SMB distribution
  SMB = SMBInitVal(SMBInitPerm(nodenumber))* exp(DeltaAcc &
         *(TSurfVal(TSurfPerm(nodenumber))- TInitVal(TInitPerm(nodenumber))))
   !PRINT *, 'SMBNEW', SMB
  !TSurf = TInitVal(TInitPerm(nodenumber)) - LapseRate * DeltaSurf + DeltaTSurf
END FUNCTION GetSMB

