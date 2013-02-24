/*
 VFHeat_SNESFEM.c
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFHeat.h"
#include "VFWell.h"
#include "VFHeat_SNESFEM.h"

/*
 VFHeat_SNESFEM
 */

/*
 ################################################################################################################
 SNES ROUTINE
 ################################################################################################################
 */
#undef __FUNCT__
#define __FUNCT__ "FEMSNESHeatSolverInitialize"
extern PetscErrorCode FEMSNESHeatSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->heatsolver = HEATSOLVER_SNESFEM;
		ctx->Hunits    = UnitaryUnits;
		ierr          = PetscOptionsEnum("-heatunits","\n\tHEAT solver","",FlowUnitName,(PetscEnum)ctx->Hunits,(PetscEnum*)&ctx->Hunits,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);	
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->KT);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->KTlhs);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATSEQAIJ,&ctx->JacT);CHKERRQ(ierr);
	} else {
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->KT);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->KTlhs);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daScal,MATMPIAIJ,&ctx->JacT);CHKERRQ(ierr);
	}
	ierr = MatZeroEntries(ctx->KT);CHKERRQ(ierr);
	ierr = MatZeroEntries(ctx->KTlhs);CHKERRQ(ierr);
	ierr = MatZeroEntries(ctx->JacT);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KTlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->JacT,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PreHeatFields);CHKERRQ(ierr);
	ierr = VecSet(ctx->PreHeatFields,0.0);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHST);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->HeatFunct);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->RHST,"RHS vector of heat equation");CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->HeatFunct,"RHS of SNES heat solver");CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSTpre);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSTpre,0.);CHKERRQ(ierr);
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesT);CHKERRQ(ierr);
	ierr = BCTInit(&ctx->bcT[0],ctx);
	ierr = BCqInit(&ctx->bcq[0],ctx);	
//	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
//	ierr = ReSETFlowBC(&ctx->bcP[0],&ctx->bcQ[0],ctx->flowcase);CHKERRQ(ierr);	
//	ierr = ReSETSourceTerms(ctx->Source,ctx->flowprop);
//	ierr = ReSETBoundaryTerms(ctx,fields);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FEMSNESHeatSolverFinalize"
extern PetscErrorCode FEMSNESHeatSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = MatDestroy(&ctx->KT);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KTlhs);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->JacT);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PreHeatFields);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHST);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->HeatFunct);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSTpre);CHKERRQ(ierr);
	ierr = SNESDestroy(&ctx->snesT);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "HeatFEMSNESSolve"
extern PetscErrorCode HeatFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	KSPConvergedReason reason;	
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           its;
	
	PetscFunctionBegin;	

	PetscFunctionReturn(0);
}





