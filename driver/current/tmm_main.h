#define MAXNUMTRACERS 100

extern PetscScalar deltaTClock, time0;
extern PetscInt maxSteps, Iter0;
extern StepTimer writeTimer;
extern PetscInt *gIndices, gLow, gHigh;
extern PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
extern PetscBool rescaleForcing;

// from here

extern PetscInt numTracers, n;
extern Vec templateVec;
extern Vec *v, *vtmp;

/* TM's */
extern Mat Ae, Ai;
extern PeriodicMat Aep, Aip;
extern TimeDependentMat Aetd, Aitd;
extern char mateFile[PETSC_MAX_PATH_LEN], matiFile[PETSC_MAX_PATH_LEN], rfsFile[PETSC_MAX_PATH_LEN];
extern PetscBool periodicMatrix;
extern PetscBool timeDependentMatrix;
extern PeriodicTimer matrixPeriodicTimer;
extern TimeDependentTimer matrixTimeDependentTimer;

/* Forcing */
extern Vec *uf, *uef;
extern PeriodicVec up[MAXNUMTRACERS];
extern char *forcingFile[MAXNUMTRACERS];
extern PeriodicTimer forcingTimer;
extern PetscInt numForcing;
extern Vec **utdf;
extern PetscScalar *tdfT;         /* array for time dependent (nonperiodic) forcing */
extern PetscScalar tf0, tf1;
extern PetscInt forcingFromFileCutOffStep;
extern PetscInt externalForcingCutOffStep;

/* Rescale forcing */
extern Vec Rfs;
extern PeriodicVec Rfsp;

/* BC's */
extern Vec *bcc, *bcf;
extern PeriodicVec bcp[MAXNUMTRACERS];
extern char *bcFile[MAXNUMTRACERS];
extern PetscInt numBC;
extern Vec **bctd;
extern PetscScalar *tdbcT;         /* array for time dependent (nonperiodic) forcing */
extern PeriodicTimer bcTimer;
extern PetscScalar tbc0, tbc1;
extern PetscInt bcCutOffStep;
extern Mat Be, Bi;
extern PeriodicMat Bep, Bip;
extern TimeDependentMat Betd, Bitd;
extern char matbeFile[PETSC_MAX_PATH_LEN], matbiFile[PETSC_MAX_PATH_LEN];
extern Vec bcTemplateVec;

/* I/O   */
extern char *iniFile[MAXNUMTRACERS];
extern char *outFile[MAXNUMTRACERS];
extern char *bcoutFile[MAXNUMTRACERS];
extern char *ufoutFile[MAXNUMTRACERS], *uefoutFile[MAXNUMTRACERS];
extern char pickupFile[PETSC_MAX_PATH_LEN];
extern char pickupoutFile[PETSC_MAX_PATH_LEN];
extern PetscBool writePickup;
extern StepTimer pickupTimer;
extern char outTimeFile[PETSC_MAX_PATH_LEN];
extern PetscBool appendOutput;
extern PetscFileMode OUTPUT_FILE_MODE;
extern PetscBool doWriteBC;
extern PetscBool doWriteUF;
extern PetscBool doWriteUEF;
extern PetscBool pickupFromFile;
extern PetscBool doTimeAverage;
extern StepTimer avgTimer;
extern char *avgOutFile[MAXNUMTRACERS];
extern Vec *vavg;
extern PetscViewer fdavgout[MAXNUMTRACERS];
extern PetscBool avgAppendOutput;
extern PetscFileMode AVG_FILE_MODE;
extern FILE *avgfptime;
extern char avgOutTimeFile[PETSC_MAX_PATH_LEN];
extern char *bcavgOutFile[MAXNUMTRACERS];
extern Vec *bcavg;
extern PetscViewer fdbcavgout[MAXNUMTRACERS];
extern char *ufavgOutFile[MAXNUMTRACERS], *uefavgOutFile[MAXNUMTRACERS];
extern Vec *ufavg, *uefavg;
extern PetscViewer fdufavgout[MAXNUMTRACERS], fduefavgout[MAXNUMTRACERS];
extern FILE *fptime;
extern PetscViewer fd, fdp, fdout[MAXNUMTRACERS];
extern PetscViewer fdbcout[MAXNUMTRACERS], fdufout[MAXNUMTRACERS], fduefout[MAXNUMTRACERS];

 #if defined (FORSPINUP) || defined (FORJACOBIAN)
extern PetscViewer fdin[MAXNUMTRACERS];
extern PetscInt itjac;
extern int fp;
 #endif
/* run time options */
extern PetscBool useExternalForcing;
extern PetscBool useForcingFromFile;
extern PetscBool usePrescribedBC;
extern PetscBool applyExternalForcing;
extern PetscBool applyForcingFromFile;
extern PetscBool applyBC;
extern PetscBool periodicForcing;
extern PetscBool timeDependentForcing;
extern PetscBool constantForcing;
extern PetscBool periodicBC;
extern PetscBool timeDependentBC;
extern PetscBool constantBC;
extern PetscBool doCalcBC;
extern PetscBool useMonitor;

extern PetscMPIInt numProcessors, myId;
extern PetscErrorCode ierr;
extern PetscBool flg1,flg2;
extern PetscScalar t1, t2, tc, tf;
extern PetscInt iLoop, Iterc;
extern PetscInt it;
extern PetscInt itr, maxValsToRead;
extern char tmpFile[PETSC_MAX_PATH_LEN];
extern PetscScalar zero, one;
extern PetscInt il;
