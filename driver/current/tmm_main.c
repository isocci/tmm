/*
   Include "petscmat.h" so that we can use matrices.
   automatically includes:
     petsc.h       - base PETSc routines   petscvec.h    - vectors
     petscsys.h    - system routines       petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"
#include "tmm_external_forcing.h"
#include "tmm_external_bc.h"
#include "tmm_timer.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_monitor.h"
#include "tmm_misfit.h"
#include "tmm_main.h"
#include "tmm_initialize.h"

PetscScalar deltaTClock, time0;
PetscInt maxSteps, Iter0;
StepTimer writeTimer;
PetscInt *gIndices, gLow, gHigh;
PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
PetscBool doMisfit = PETSC_FALSE;
PetscBool rescaleForcing = PETSC_FALSE;

extern PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers,
                                  PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool usePrescribedBC,
                                  Vec *v, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp);

extern PetscErrorCode iniTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v, PetscBool append);
extern PetscErrorCode doTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode finalizeTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr);

extern PetscErrorCode waitForSignal(PetscInt waitTime);

#undef __FUNCT__
#define __FUNCT__ "main"


int main(int argc,char **args)
{
        //----------------------------------------------------------------------------

        PetscInt numTracers, n;
        Vec templateVec;
        Vec *v, *vtmp;
/* TM's */
        Mat Ae, Ai;
//  PeriodicMat Aep, Aip;
//  TimeDependentMat Aetd, Aitd;
        char mateFile[PETSC_MAX_PATH_LEN], matiFile[PETSC_MAX_PATH_LEN], rfsFile[PETSC_MAX_PATH_LEN];
        PetscBool periodicMatrix = PETSC_FALSE;
        PetscBool timeDependentMatrix = PETSC_FALSE;
        PeriodicTimer matrixPeriodicTimer;
        TimeDependentTimer matrixTimeDependentTimer;

/* Forcing */
        Vec *uf, *uef;
        PeriodicVec up[MAXNUMTRACERS];
        char *forcingFile[MAXNUMTRACERS];
        PeriodicTimer forcingTimer;
        PetscInt numForcing;
        Vec **utdf;
        PetscScalar *tdfT; /* array for time dependent (nonperiodic) forcing */
        PetscScalar tf0, tf1;
        PetscInt forcingFromFileCutOffStep = -1;
        PetscInt externalForcingCutOffStep = -1;

/* Rescale forcing */
        Vec Rfs;
        PeriodicVec Rfsp;

/* BC's */
        Vec *bcc, *bcf;
        PeriodicVec bcp[MAXNUMTRACERS];
        char *bcFile[MAXNUMTRACERS];
        PetscInt numBC;
        Vec **bctd;
        PetscScalar *tdbcT; /* array for time dependent (nonperiodic) forcing */
        PeriodicTimer bcTimer;
        PetscScalar tbc0, tbc1;
        PetscInt bcCutOffStep = -1;
        Mat Be, Bi;
        PeriodicMat Bep, Bip;
        TimeDependentMat Betd, Bitd;
        char matbeFile[PETSC_MAX_PATH_LEN], matbiFile[PETSC_MAX_PATH_LEN];
        Vec bcTemplateVec;

/* I/O   */
        char *iniFile[MAXNUMTRACERS];
        char *outFile[MAXNUMTRACERS];
//  char *bcoutFile[MAXNUMTRACERS];
        char *ufoutFile[MAXNUMTRACERS]; // *uefoutFile[MAXNUMTRACERS];
        char pickupFile[PETSC_MAX_PATH_LEN];
        char pickupoutFile[PETSC_MAX_PATH_LEN];
        PetscBool writePickup = PETSC_FALSE;
        StepTimer pickupTimer;
        char outTimeFile[PETSC_MAX_PATH_LEN];
        PetscBool appendOutput = PETSC_FALSE;
        PetscFileMode OUTPUT_FILE_MODE;
        PetscBool doWriteBC = PETSC_FALSE;
        PetscBool doWriteUF = PETSC_FALSE;
        PetscBool doWriteUEF = PETSC_FALSE;
        PetscBool pickupFromFile = PETSC_FALSE;
        PetscBool doTimeAverage = PETSC_FALSE;
        StepTimer avgTimer;
        char *avgOutFile[MAXNUMTRACERS];
        Vec *vavg;
        PetscViewer fdavgout[MAXNUMTRACERS];
        PetscBool avgAppendOutput = PETSC_FALSE;
        PetscFileMode AVG_FILE_MODE;
        FILE *avgfptime;
        char avgOutTimeFile[PETSC_MAX_PATH_LEN];
        char *bcavgOutFile[MAXNUMTRACERS];
        Vec *bcavg;
        PetscViewer fdbcavgout[MAXNUMTRACERS];
        char *ufavgOutFile[MAXNUMTRACERS], *uefavgOutFile[MAXNUMTRACERS];
        Vec *ufavg, *uefavg;
        PetscViewer fdufavgout[MAXNUMTRACERS], fduefavgout[MAXNUMTRACERS];
        FILE *fptime;
        PetscViewer fd, fdp, fdout[MAXNUMTRACERS];
        PetscViewer fdbcout[MAXNUMTRACERS], fdufout[MAXNUMTRACERS], fduefout[MAXNUMTRACERS];

#if defined (FORSPINUP) || defined (FORJACOBIAN)
        PetscViewer fdin[MAXNUMTRACERS];
        PetscInt itjac;
        int fp;
#endif

/* run time options */
        PetscBool useExternalForcing = PETSC_FALSE;
        PetscBool useForcingFromFile = PETSC_FALSE;
        PetscBool usePrescribedBC = PETSC_FALSE;
//  PetscBool applyExternalForcing = PETSC_FALSE;
//  PetscBool applyForcingFromFile = PETSC_FALSE;
        PetscBool applyBC = PETSC_FALSE;
        PetscBool periodicForcing = PETSC_FALSE;
        PetscBool timeDependentForcing = PETSC_FALSE;
//  PetscBool constantForcing = PETSC_FALSE;
        PetscBool periodicBC = PETSC_FALSE;
        PetscBool timeDependentBC = PETSC_FALSE;
        PetscBool constantBC = PETSC_FALSE;
        PetscBool doCalcBC = PETSC_FALSE;
        PetscBool useMonitor = PETSC_FALSE;

        PetscMPIInt numProcessors, myId;
        PetscErrorCode ierr;
        PetscBool flg1,flg2;
        PetscScalar t1, t2, tc, tf;
        PetscInt iLoop, Iterc;
        PetscInt it;
        PetscInt itr; // maxValsToRead;
        char tmpFile[PETSC_MAX_PATH_LEN];
        PetscScalar zero = 0.0, one = 1.0;
        PetscInt il;
        //----------------------------------------------------------------------------

        int status = initialize(argc, **args);

/* Start time stepping loop */
        ierr = PetscTime(&t1); CHKERRQ(ierr); /* start counting wall clock time */
        for (iLoop = 1; iLoop <= maxSteps; iLoop++) {
/*  iLoop -> iLoop+1 (convention) */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
                tc=time0 + deltaTClock*(iLoop-1); /*  current time (time at beginning of step) */
                tf=time0 + deltaTClock*iLoop; /*  future time (time at end of step) */
                Iterc=Iter0+iLoop-1;
#else
                tc=time0 + deltaTClock*(itjac-1); /*  current time (time at beginning of step) */
                tf=time0 + deltaTClock*itjac; /*  future time (time at end of step) */
                Iterc=Iter0+itjac-1;
                itjac++;
#endif

#ifdef FORJACOBIAN
                if (!doTimeAverage) {
                        for (itr=0; itr<numTracers; itr++) {
/* reread initial conditions (these were read once already above, but do it again since its easier to code) */
                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading initial condition from file %s\n", itr,iniFile[itr]); CHKERRQ(ierr);
                                ierr = VecLoad(v[itr],fdin[itr]); CHKERRQ(ierr); /* IntoVector */
                        }
                        if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef); CHKERRQ(ierr);
                }
#endif

#ifndef FORJACOBIAN
/*  interpolate Ae,Ai,uf,uef,bcc to current time (tc) and bcf to future time (tf) */
                if (periodicMatrix) {
                        ierr = interpPeriodicMatrix(tc,&Ae,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                                    matrixPeriodicTimer.tdp,&Aep,mateFile);
                        ierr = interpPeriodicMatrix(tc,&Ai,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                                    matrixPeriodicTimer.tdp,&Aip,matiFile);
                        if (rescaleForcing) {
                                ierr = interpPeriodicVector(tc,&Rfs,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                                            matrixPeriodicTimer.tdp,&Rfsp,rfsFile);
                        }
                } else if (timeDependentMatrix) {
                        ierr = interpTimeDependentMatrix(tc,&Ae,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Aetd,mateFile);
                        ierr = interpTimeDependentMatrix(tc,&Ai,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Aitd,matiFile);
                }

/*  Forcing     */
                if (applyForcingFromFile) {
                        if ((forcingFromFileCutOffStep>0) && (iLoop==(forcingFromFileCutOffStep+1))) {
                                applyForcingFromFile = PETSC_FALSE;
                                for (itr=0; itr<numTracers; itr++) {
                                        VecSet(uf[itr],zero);
                                }
                        } else {
                                if (periodicForcing) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = interpPeriodicVector(tc,&uf[itr],forcingTimer.cyclePeriod,forcingTimer.numPerPeriod,forcingTimer.tdp,&up[itr],forcingFile[itr]);
                                        }
                                } else if (timeDependentForcing) {
                                        ierr = interpTimeDependentVector(tc,uf,numTracers,numForcing,tdfT,utdf); CHKERRQ(ierr);
                                }
                                if (rescaleForcing) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecPointwiseMult(uf[itr],Rfs,uf[itr]); CHKERRQ(ierr);
                                        }
                                }
                        }
                }
#endif

                if (applyExternalForcing) {
                        if ((externalForcingCutOffStep>0) && (iLoop==(externalForcingCutOffStep+1))) {
                                applyExternalForcing = PETSC_FALSE;
                                for (itr=0; itr<numTracers; itr++) {
                                        VecSet(uef[itr],zero);
                                }
                        } else {
                                ierr = calcExternalForcing(tc,Iterc,iLoop,numTracers,v,uef); CHKERRQ(ierr); /* Compute external forcing in uef */
                                if (rescaleForcing) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecPointwiseMult(uef[itr],Rfs,uef[itr]); CHKERRQ(ierr);
                                        }
                                }
                        }
                }

#ifndef FORJACOBIAN
                if (applyBC) {
                        if ((bcCutOffStep>0) && (iLoop==(bcCutOffStep+1))) {
                                applyBC = PETSC_FALSE;
                                doCalcBC = PETSC_FALSE;
                                for (itr=0; itr<numTracers; itr++) {
                                        VecSet(bcc[itr],zero);
                                        VecSet(bcf[itr],zero);
                                }
                        } else {
                                if (periodicMatrix) {
                                        ierr = interpPeriodicMatrix(tc,&Be,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                                                    matrixPeriodicTimer.tdp,&Bep,matbeFile);
                                        ierr = interpPeriodicMatrix(tc,&Bi,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                                                    matrixPeriodicTimer.tdp,&Bip,matbiFile);
                                } else if (timeDependentMatrix) {
                                        ierr = interpTimeDependentMatrix(tc,&Be,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Betd,matbeFile);
                                        ierr = interpTimeDependentMatrix(tc,&Bi,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Bitd,matbiFile);
                                }
                                if (periodicBC) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = interpPeriodicVector(tc,&bcc[itr],bcTimer.cyclePeriod,bcTimer.numPerPeriod,bcTimer.tdp,&bcp[itr],bcFile[itr]);
                                                ierr = interpPeriodicVector(tf,&bcf[itr],bcTimer.cyclePeriod,bcTimer.numPerPeriod,bcTimer.tdp,&bcp[itr],bcFile[itr]);
                                        }
                                } else if (timeDependentBC) {
                                        ierr = interpTimeDependentVector(tc,bcc,numTracers,numBC,tdbcT,bctd); CHKERRQ(ierr);
                                        ierr = interpTimeDependentVector(tf,bcf,numTracers,numBC,tdbcT,bctd); CHKERRQ(ierr);
                                } else if (doCalcBC) {
                                        ierr = calcBC(tc,Iterc,tc+deltaTClock,Iterc+1,iLoop,numTracers,v,bcc,bcf); CHKERRQ(ierr); /* Compute BC in bcc and bcf */
                                }
                        }

                }

/*  Write out BC at first time step */
                if (doWriteBC) {
                        if ((iLoop == 1) && (!appendOutput)) {
                                for (itr=0; itr<numTracers; itr++) {
                                        ierr = VecView(bcc[itr],fdbcout[itr]); CHKERRQ(ierr);
                                }
                        }
                }

                ierr = forwardStep(tc,Iterc,deltaTClock,numTracers,applyForcingFromFile,applyExternalForcing,applyBC,
                                   v,Ae,Ai,Be,Bi,uf,uef,bcc,bcf,vtmp); CHKERRQ(ierr);
#endif
                tc=time0 + deltaTClock*iLoop; /*  time at end of step */

                if (useMonitor) {
                        ierr = updateMonitor(tc,iLoop,numTracers,v); CHKERRQ(ierr);
                        ierr = writeMonitor(tc,iLoop,numTracers,v); CHKERRQ(ierr);
                }

                if (doMisfit) {
                        ierr = calcMisfit(tc,iLoop,numTracers,v); CHKERRQ(ierr);
                        ierr = writeMisfit(tc,iLoop,numTracers,v); CHKERRQ(ierr);
                }

/* write output */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
                if (applyExternalForcing) {
                        ierr = writeExternalForcing(tc,iLoop,numTracers,v,uef); CHKERRQ(ierr);
                }
                if (doCalcBC) {
                        ierr = writeBC(tc,iLoop,numTracers,v,bcc,bcf); CHKERRQ(ierr);
                }
                if (Iter0+iLoop>=writeTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
                        if (writeTimer.count<=writeTimer.numTimeSteps) {
                                writeTimer.count++;
                        }
                        if (writeTimer.count==writeTimer.numTimeSteps) { /* time to write out */
                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0+iLoop,tc); CHKERRQ(ierr);
                                for (itr=0; itr<numTracers; itr++) {
                                        ierr = VecView(v[itr],fdout[itr]); CHKERRQ(ierr);
                                }

                                if (doWriteUF) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecView(uf[itr],fdufout[itr]); CHKERRQ(ierr);
                                        }
                                }

                                if (doWriteUEF) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecView(uef[itr],fduefout[itr]); CHKERRQ(ierr);
                                        }
                                }

                                if (doWriteBC) {
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecView(bcf[itr],fdbcout[itr]); CHKERRQ(ierr);
                                        }
                                }

                                ierr = updateStepTimer("write_", Iter0+iLoop, &writeTimer); CHKERRQ(ierr);
                        }
                }

                ierr = doTMMWrite(tc,iLoop,numTracers,v); CHKERRQ(ierr);

                if (writePickup) {
                        if (Iter0+iLoop>=pickupTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
                                if (pickupTimer.count<=pickupTimer.numTimeSteps) {
                                        pickupTimer.count++;
                                }
                                if (pickupTimer.count==pickupTimer.numTimeSteps) { /* time to write out */
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing pickup at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                        strcpy(tmpFile,"");
                                        sprintf(tmpFile,"pickup_%d.petsc",Iter0+iLoop);
                                        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fdp); CHKERRQ(ierr);
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecView(v[itr],fdp); CHKERRQ(ierr);
                                        }
                                        ierr = PetscViewerDestroy(&fdp); CHKERRQ(ierr);

                                        ierr = updateStepTimer("pickup_", Iter0+iLoop, &pickupTimer); CHKERRQ(ierr);
                                }
                        }
                }

#else
#ifdef FORSPINUP
                if (Iter0+iLoop>=writeTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
                        if (writeTimer.count<=writeTimer.numTimeSteps) {
                                writeTimer.count++;
                        }
                        if (writeTimer.count==writeTimer.numTimeSteps) { /* time to write out */
                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                for (itr=0; itr<numTracers; itr++) {
                                        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]); CHKERRQ(ierr);
                                        ierr = VecView(v[itr],fdout[itr]); CHKERRQ(ierr);
                                        ierr = PetscViewerDestroy(&fdout[itr]); CHKERRQ(ierr);
                                }

                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n"); CHKERRQ(ierr);
                                ierr = waitForSignal(10); CHKERRQ(ierr);

                                for (itr=0; itr<numTracers; itr++) {
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]); CHKERRQ(ierr);
                                        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                        ierr = VecLoad(v[itr],fd); CHKERRQ(ierr); /* IntoVector */
                                        ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
                                }

                                ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                ierr = PetscViewerBinaryGetDescriptor(fd,&fp); CHKERRQ(ierr);
                                ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT); CHKERRQ(ierr);
                                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac); CHKERRQ(ierr);

                                tc=time0 + deltaTClock*(itjac-1); /*  current time (time at beginning of step) */
                                tf=time0 + deltaTClock*itjac; /*  future time (time at end of step) */
                                Iterc=Iter0+itjac-1;

                                if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef); CHKERRQ(ierr);
                                if (doCalcBC) ierr = reInitializeCalcBC(tc,Iterc,tc+deltaTClock,Iterc+1,numTracers,v,bcc,bcf); CHKERRQ(ierr);
                        }
                }
#endif
#ifdef FORJACOBIAN
                if (doTimeAverage) {
                        if (Iter0+iLoop>=avgTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
                                if (avgTimer.count<=avgTimer.numTimeSteps) { /* still within same averaging block so accumulate */
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecAXPY(vavg[itr],one,uef[itr]); CHKERRQ(ierr);
                                        }
                                        avgTimer.count++;
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating: %d\n", iLoop);CHKERRQ(ierr);                       */
                                }
                                if (avgTimer.count==avgTimer.numTimeSteps) { /* time to write averages to file */
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average q at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecScale(vavg[itr],1.0/avgTimer.count); CHKERRQ(ierr);
                                                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]); CHKERRQ(ierr);
                                                ierr = VecView(vavg[itr],fdout[itr]); CHKERRQ(ierr);
                                                ierr = PetscViewerDestroy(&fdout[itr]); CHKERRQ(ierr);
                                                ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
                                        }
                                        ierr = updateStepTimer("avg_", Iter0+iLoop, &avgTimer); CHKERRQ(ierr);

                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n"); CHKERRQ(ierr);
                                        ierr = waitForSignal(10); CHKERRQ(ierr);

                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]); CHKERRQ(ierr);
                                                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                                ierr = VecLoad(v[itr],fd); CHKERRQ(ierr); /* IntoVector */
                                                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
                                        }

                                        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                        ierr = PetscViewerBinaryGetDescriptor(fd,&fp); CHKERRQ(ierr);
                                        ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT); CHKERRQ(ierr);
                                        ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac); CHKERRQ(ierr);

                                        tc=time0 + deltaTClock*(itjac-1); /*  current time (time at beginning of step) */
                                        tf=time0 + deltaTClock*itjac; /*  future time (time at end of step) */
                                        Iterc=Iter0+itjac-1;

                                        if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef); CHKERRQ(ierr);

                                } else {

                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]); CHKERRQ(ierr);
                                                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                                ierr = VecLoad(v[itr],fd); CHKERRQ(ierr); /* IntoVector */
                                                ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
                                        }

                                        tc=time0 + deltaTClock*(itjac-1); /* this is now the time at end of time step */
                                        Iterc=Iter0+itjac-1;
                                        if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef); CHKERRQ(ierr);
                                }
                        }
                } else { /* no time averaging */
                        if (Iter0+iLoop>=writeTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
                                if (writeTimer.count<=writeTimer.numTimeSteps) {
                                        writeTimer.count++;
                                }
                                if (writeTimer.count==writeTimer.numTimeSteps) { /* time to write out */
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing q at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecView(uef[itr],fdout[itr]); CHKERRQ(ierr);
                                        }
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = PetscViewerDestroy(&fdout[itr]); CHKERRQ(ierr);
                                                ierr = PetscViewerDestroy(&fdin[itr]); CHKERRQ(ierr);
                                        }
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n"); CHKERRQ(ierr);
                                        ierr = waitForSignal(10); CHKERRQ(ierr);
                                        /*      open files here for I/O */
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]); CHKERRQ(ierr);
                                                ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fdin[itr]); CHKERRQ(ierr);
                                        }

                                        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd); CHKERRQ(ierr);
                                        ierr = PetscViewerBinaryGetDescriptor(fd,&fp); CHKERRQ(ierr);
                                        ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT); CHKERRQ(ierr);
                                        ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac); CHKERRQ(ierr);
                                }
                        }
                }
#endif
#endif
#ifndef FORJACOBIAN
                if (doTimeAverage) {
                        if (Iter0+iLoop>=avgTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
                                if (avgTimer.count<=avgTimer.numTimeSteps) { /* still within same averaging block so accumulate */
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating for time average\n");CHKERRQ(ierr);               */
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecAXPY(vavg[itr],one,v[itr]); CHKERRQ(ierr);
                                        }
                                        if (doWriteUF) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecAXPY(ufavg[itr],one,uf[itr]); CHKERRQ(ierr);
                                                }
                                        }
                                        if (doWriteUEF) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecAXPY(uefavg[itr],one,uef[itr]); CHKERRQ(ierr);
                                                }
                                        }
                                        if (doWriteBC) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecAXPY(bcavg[itr],one,bcf[itr]); CHKERRQ(ierr);
                                                }
                                        }
                                        avgTimer.count++;
                                }
                                if (avgTimer.count==avgTimer.numTimeSteps) { /* time to write averages to file */
                                        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average at time %10.5f, step %d\n", tc, Iter0+iLoop); CHKERRQ(ierr);
                                        ierr = PetscFPrintf(PETSC_COMM_WORLD,avgfptime,"%d   %10.5f\n",Iter0+iLoop,tc); CHKERRQ(ierr);
                                        for (itr=0; itr<numTracers; itr++) {
                                                ierr = VecScale(vavg[itr],1.0/avgTimer.count); CHKERRQ(ierr);
                                                ierr = VecView(vavg[itr],fdavgout[itr]); CHKERRQ(ierr);
                                                ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
                                        }
                                        if (doWriteUF) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecScale(ufavg[itr],1.0/avgTimer.count); CHKERRQ(ierr);
                                                        ierr = VecView(ufavg[itr],fdufavgout[itr]); CHKERRQ(ierr);
                                                        ierr = VecSet(ufavg[itr],zero); CHKERRQ(ierr);
                                                }
                                        }
                                        if (doWriteUEF) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecScale(uefavg[itr],1.0/avgTimer.count); CHKERRQ(ierr);
                                                        ierr = VecView(uefavg[itr],fduefavgout[itr]); CHKERRQ(ierr);
                                                        ierr = VecSet(uefavg[itr],zero); CHKERRQ(ierr);
                                                }
                                        }
                                        if (doWriteBC) {
                                                for (itr=0; itr<numTracers; itr++) {
                                                        ierr = VecScale(bcavg[itr],1.0/avgTimer.count); CHKERRQ(ierr);
                                                        ierr = VecView(bcavg[itr],fdbcavgout[itr]); CHKERRQ(ierr);
                                                        ierr = VecSet(bcavg[itr],zero); CHKERRQ(ierr);
                                                }
                                        }
                                        ierr = updateStepTimer("avg_", Iter0+iLoop, &avgTimer); CHKERRQ(ierr);
                                }
                        }
                }
#endif
        } /* end of time-stepping loop */

#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
        for (itr=0; itr<numTracers; itr++) {
                ierr = PetscViewerDestroy(&fdout[itr]); CHKERRQ(ierr);
        }
        ierr = PetscFClose(PETSC_COMM_WORLD,fptime); CHKERRQ(ierr);

        ierr = finalizeTMMWrite(tc,maxSteps,numTracers); CHKERRQ(ierr);

        if (doWriteUF) {
                for (itr=0; itr<numTracers; itr++) {
                        ierr = PetscViewerDestroy(&fdufout[itr]); CHKERRQ(ierr);
                }
        }

        if (doWriteUEF) {
                for (itr=0; itr<numTracers; itr++) {
                        ierr = PetscViewerDestroy(&fduefout[itr]); CHKERRQ(ierr);
                }
        }

        if (doWriteBC) {
                for (itr=0; itr<numTracers; itr++) {
                        ierr = PetscViewerDestroy(&fdbcout[itr]); CHKERRQ(ierr);
                }
        }
#endif

        if (doTimeAverage) {
                ierr = PetscFClose(PETSC_COMM_WORLD,avgfptime); CHKERRQ(ierr);
                for (itr=0; itr<numTracers; itr++) {
                        ierr = PetscViewerDestroy(&fdavgout[itr]); CHKERRQ(ierr);
                }
                ierr = VecDestroyVecs(numTracers,&vavg); CHKERRQ(ierr);
                if (doWriteUF) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = PetscViewerDestroy(&fdufavgout[itr]); CHKERRQ(ierr);
                        }
                        ierr = VecDestroyVecs(numTracers,&ufavg); CHKERRQ(ierr);
                }
                if (doWriteUEF) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = PetscViewerDestroy(&fduefavgout[itr]); CHKERRQ(ierr);
                        }
                        ierr = VecDestroyVecs(numTracers,&uefavg); CHKERRQ(ierr);
                }
                if (doWriteBC) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = PetscViewerDestroy(&fdbcavgout[itr]); CHKERRQ(ierr);
                        }
                        ierr = VecDestroyVecs(numTracers,&bcavg); CHKERRQ(ierr);
                }
        }

/* write final pickup */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,pickupoutFile,FILE_MODE_WRITE,&fdp); CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
                ierr = VecView(v[itr],fdp); CHKERRQ(ierr);
        }
        ierr = PetscViewerDestroy(&fdp); CHKERRQ(ierr);

        ierr=PetscTime(&t2); CHKERRQ(ierr); /* stop counting wall clock time */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Wall clock time: %10.5f\n", t2-t1); CHKERRQ(ierr);

        /* Free data structures */
        ierr = VecDestroyVecs(numTracers,&v); CHKERRQ(ierr);
        ierr = VecDestroyVecs(numTracers,&vtmp); CHKERRQ(ierr);
        ierr = MatDestroy(&Ae); CHKERRQ(ierr);
        ierr = MatDestroy(&Ai); CHKERRQ(ierr);

        if (periodicMatrix) {
                ierr = destroyPeriodicMat(&Aep); CHKERRQ(ierr);
                ierr = destroyPeriodicMat(&Aip); CHKERRQ(ierr);
        } else if (timeDependentMatrix) {
                ierr = destroyTimeDependentMat(&Aetd);
                ierr = destroyTimeDependentMat(&Aitd);
        }

        if (rescaleForcing) {
                ierr = VecDestroy(&Rfs); CHKERRQ(ierr);
                if (periodicMatrix) {
                        ierr = destroyPeriodicVec(&Rfsp); CHKERRQ(ierr);
                }
        }

        if (useMonitor) {
                ierr = finalizeMonitor(tc,maxSteps,numTracers); CHKERRQ(ierr);
        }

        if (doMisfit) {
                ierr = finalizeMisfit(tc,maxSteps,numTracers); CHKERRQ(ierr);
        }

        if (useExternalForcing) {
                ierr = VecDestroyVecs(numTracers,&uef); CHKERRQ(ierr);
                ierr = finalizeExternalForcing(tc,maxSteps,numTracers); CHKERRQ(ierr);
        }

        if (useForcingFromFile) {
                ierr = VecDestroyVecs(numTracers,&uf); CHKERRQ(ierr);
                if (periodicForcing) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = destroyPeriodicVec(&up[itr]); CHKERRQ(ierr);
                        }
                } else if (timeDependentForcing) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = VecDestroyVecs(numForcing,&utdf[itr]); CHKERRQ(ierr);
                        }
                }
        }

        if (usePrescribedBC) {
                ierr = MatDestroy(&Be); CHKERRQ(ierr);
                ierr = MatDestroy(&Bi); CHKERRQ(ierr);
                if (periodicMatrix) {
                        ierr = destroyPeriodicMat(&Bep); CHKERRQ(ierr);
                        ierr = destroyPeriodicMat(&Bip); CHKERRQ(ierr);
                } else if (timeDependentMatrix) {
                        ierr = destroyTimeDependentMat(&Betd);
                        ierr = destroyTimeDependentMat(&Bitd);
                }
                ierr = VecDestroyVecs(numTracers,&bcc); CHKERRQ(ierr);
                ierr = VecDestroyVecs(numTracers,&bcf); CHKERRQ(ierr);
                if (periodicBC) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = destroyPeriodicVec(&bcp[itr]); CHKERRQ(ierr);
                        }
                } else if (timeDependentBC) {
                        for (itr=0; itr<numTracers; itr++) {
                                ierr = VecDestroyVecs(numBC,&bctd[itr]); CHKERRQ(ierr);
                        }
                } else if (doCalcBC) {
                        ierr = finalizeCalcBC(tc,maxSteps,numTracers); CHKERRQ(ierr);
                        ierr = PetscFree(gBCIndices); CHKERRQ(ierr);
                }
        }

        ierr = PetscFree(gIndices); CHKERRQ(ierr);

        ierr = PetscFinalize(); CHKERRQ(ierr);
        return 0;
}
