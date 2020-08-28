extern int initialize(int, char **);

extern PetscErrorCode iniTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v, PetscBool append);

//global in main---------------------------------------------------------------
extern PetscScalar deltaTClock, time0;
extern PetscInt maxSteps, Iter0;
extern StepTimer writeTimer;
extern PetscInt *gIndices, gLow, gHigh;
extern PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
extern PetscBool doMisfit;
extern PetscBool rescaleForcing;
//------------------------------------------------------------------------------

//variables made global in initialize-------------------------------------------
//extern Mat Ae, Ai;
extern PeriodicMat Aep, Aip;
extern TimeDependentMat Aetd, Aitd;

extern Vec **utdf;
extern PetscInt numForcing;
extern char *forcingFile[MAXNUMTRACERS];

extern PeriodicVec bcp[MAXNUMTRACERS];
extern Vec **bctd;
extern Mat Be, Bi;
extern PeriodicMat Bep, Bip;
extern TimeDependentMat Betd, Bitd;

extern PeriodicVec up[MAXNUMTRACERS];

extern PetscBool applyForcingFromFile;
extern PetscBool applyExternalForcing;

extern PeriodicVec Rfsp;

extern PetscInt numBC;
//------------------------------------------------------------------------------
