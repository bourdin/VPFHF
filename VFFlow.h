/*
   VFFlow.h
   Generic interface to flow solvers

     (c) 2010-2012 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#ifndef VFFLOW_H
#define VFFLOW_H


typedef struct {
	PetscReal       mu;			/* Fluid viscosity						*/
	PetscReal       rho;        /* Fluid density						*/
	PetscReal       cf;         /* Fluid compressibility				*/
	PetscReal		betac;		/* Conversion constant					*/
	PetscReal		gammac;		/*Conversion parameter					*/
	PetscReal		alphac;		/*Conversion parameter					*/
	PetscReal		g;
} FluidProp; //change them to Vec later

extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BCPInit(BC *BCP,VFCtx *ctx);
extern PetscErrorCode BCQInit(BC *BCQ,VFCtx *ctx);
extern PetscErrorCode VecApplyPressureBC_FEM(Vec RHS,Vec BCF,BC *BC);

#endif /* VFFLOW_H */

