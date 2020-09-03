#ifndef EXTERN
#define EXTERN extern
#endif


//originally global in main-----------------------------------------------------
EXTERN PetscScalar deltaTClock, time0;
EXTERN PetscInt maxSteps, Iter0;
EXTERN StepTimer writeTimer;
EXTERN PetscInt *gIndices, gLow, gHigh;
EXTERN PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
EXTERN PetscBool doMisfit;
EXTERN PetscBool rescaleForcing;
//------------------------------------------------------------------------------

//variables made global---------------------------------------------------------
//EXTERN Mat Ae, Ai;
EXTERN FILE *fptime;

EXTERN PeriodicMat Aep, Aip;
EXTERN TimeDependentMat Aetd, Aitd;

EXTERN Vec **utdf;
EXTERN PetscInt numForcing;
EXTERN char *forcingFile[MAXNUMTRACERS];

EXTERN PeriodicVec bcp[MAXNUMTRACERS];
EXTERN Vec **bctd;
EXTERN Mat Be, Bi;
EXTERN PeriodicMat Bep, Bip;
EXTERN TimeDependentMat Betd, Bitd;

EXTERN PeriodicVec up[MAXNUMTRACERS];

EXTERN PetscBool applyForcingFromFile;
EXTERN PetscBool applyExternalForcing;

EXTERN PeriodicVec Rfsp;

EXTERN PetscInt numBC;

EXTERN char *iniFile[MAXNUMTRACERS];
EXTERN char *outFile[MAXNUMTRACERS];
EXTERN char *bcoutFile[MAXNUMTRACERS];
EXTERN char *ufoutFile[MAXNUMTRACERS], *uefoutFile[MAXNUMTRACERS];

EXTERN Vec *bcc, *bcf;
EXTERN char *bcFile[MAXNUMTRACERS];
EXTERN PetscScalar *tdbcT; /* array for time dependent (nonperiodic) forcing */

EXTERN char *avgOutFile[MAXNUMTRACERS];
EXTERN Vec *vavg;

EXTERN FILE *avgfptime;
EXTERN char avgOutTimeFile[PETSC_MAX_PATH_LEN];
EXTERN char *bcavgOutFile[MAXNUMTRACERS];
EXTERN Vec *bcavg;

EXTERN char *ufavgOutFile[MAXNUMTRACERS], *uefavgOutFile[MAXNUMTRACERS];
EXTERN Vec *ufavg, *uefavg;

EXTERN Vec *v, *vtmp;

EXTERN Vec *uf, *uef;
EXTERN PetscScalar *tdfT; /* array for time dependent (nonperiodic) forcing */
