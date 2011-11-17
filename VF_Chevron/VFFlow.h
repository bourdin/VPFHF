/*
   VFFlow.c
   Generic interface to flow solvers

     (c) 2010-2011 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"

/*
 VFFlow.h
 VF - Fluid Flow
 */

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"



typedef enum {
	NORMALVELOCITY,
	PRESSURE
} FlowBCTYPE;

typedef struct {
	FlowBCTYPE     face[6];
	FlowBCTYPE     edge[12];
	FlowBCTYPE     vertex[8];   
} FLOWBC;

static const char *FLOWBCTYPE_NAME[] = {
"NORMALVELOCITY",
"PRESSURE",
"FLOWBCTYPE_NAME",
"",
0
};

typedef enum { 
	ALLNORMALFLOWBC,
	ALLPRESSUREBC
} FlowCases;

typedef struct {
	PetscReal       mu;			/* Fluid viscosity						*/
	PetscReal       rho;        /* Fluid density						*/
	PetscReal       cf;         /* Fluid compressibility				*/
	PetscReal		betac;		/* Conversion constant					*/
	PetscReal		gammac;		/*Conversion parameter					*/
	PetscReal		alphac;		/*Conversion parameter					*/
	PetscReal		g;
} FluidProp; //change them to Vec later

extern PetscErrorCode GetFlowBC(FLOWBC *BC, FlowCases setbc)
extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields)
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx)
#endif /* VFFLOW_H */



/*
 1.	Flowsolver Initialize
 2.	Get boundary conditions
 3.	Mattrix Assemble
 4.	Matrix apply boundary conditions
 5.	Vector apply boundary conditions
 6.	KSp solve
 7.	Finalize
 
 */

































#endif /* VFFLOW_H */
