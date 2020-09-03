#define MAXNUMTRACERS 100

extern PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers,
                                  PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool usePrescribedBC,
                                  Vec *v, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp);

extern PetscErrorCode iniTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v, PetscBool append);
extern PetscErrorCode doTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode finalizeTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr);

extern PetscErrorCode waitForSignal(PetscInt waitTime);
