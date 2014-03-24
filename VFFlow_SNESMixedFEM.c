/*
 VFFlow_MixedFEM.c
 A mixed finite elements Darcy solver based on the method in
 Masud, A. and Hughes, T. J. (2002). A stabilized mixed finite element method for
 Darcy flow. Computer Methods in Applied Mechanics and Engineering, 191(3940):43414370.
 
 (c) 2011-2012 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_SNESMixedFEM.h"
#include "VFFlow_KSPMixedFEM.h"
#include "VFWell.h"

/*
 VFFlow_DarcyMixedFEMSNES
 */

/*
 ################################################################################################################
 SNES ROUTINE
 ################################################################################################################
 */
#undef __FUNCT__
#define __FUNCT__ "MixedFEMSNESFlowSolverInitialize"
extern PetscErrorCode MixedFEMSNESFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
  PetscFunctionBegin;
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesVelP);CHKERRQ(ierr);
  ierr = SNESAppendOptionsPrefix(ctx->snesVelP,"FlowSnes_");CHKERRQ(ierr);
  ierr = SNESSetFromOptions(ctx->snesVelP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMSNESFlowSolverFinalize"
extern PetscErrorCode MixedFEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;
	ierr = SNESDestroy(&ctx->snesVelP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESIFunction"
extern PetscErrorCode FormSNESIFunction(SNES snes,Vec VelnPress,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			 *ctx=(VFCtx*)user;
	Vec         VecRHS;
	PetscReal		theta;
	PetscReal		one_minus_theta;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	theta = ctx->flowprop.theta;
	one_minus_theta = (1.-theta);
	ierr = VecDuplicate(ctx->RHSVelP,&VecRHS);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,VecRHS);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,one_minus_theta,theta,ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = MatMultAdd(ctx->KVelPlhs,ctx->PreFlowFields,VecRHS,VecRHS);CHKERRQ(ierr);
	ierr = VecApplyVelocityBC(VecRHS,ctx->VelBCArray,&ctx->bcQ[0],ctx);CHKERRQ(ierr);
	ierr = MatMult(ctx->KVelP,VelnPress,Func);CHKERRQ(ierr);
  ierr = VecAXPY(Func,-1.0,VecRHS);CHKERRQ(ierr);
	ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFlowFEMSNESSolve"
extern PetscErrorCode MixedFlowFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	SNESConvergedReason reason;	
	PetscReal          ****VelnPress_array;
	PetscReal          ***Press_array;
	PetscReal          ****vel_array;
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           its;
	PetscReal           Velmin,Velmax;
	PetscReal           Pmin,Pmax;

	PetscFunctionBegin;
	ierr = VecCopy(fields->V,ctx->V);CHKERRQ(ierr);
	ierr = VecCopy(fields->U,ctx->U);CHKERRQ(ierr);
  ierr = FlowMatnVecAssemble(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,fields,ctx);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesVelP,ctx->FlowFunct,FormSNESIFunction,ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(ctx->snesVelP,ctx->JacVelP,ctx->JacVelP,FormSNESIJacobian,ctx);CHKERRQ(ierr);
	if (ctx->verbose > 1) {
		ierr = SNESMonitorSet(ctx->snesVelP,FEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	}
    ierr = SNESSolve(ctx->snesVelP,PETSC_NULL,fields->VelnPress);CHKERRQ(ierr);
	ierr = SNESGetConvergedReason(ctx->snesVelP,&reason);CHKERRQ(ierr);
	if (reason < 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] snes_MixedFlowSolver diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = SNESGetIterationNumber(ctx->snesVelP,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      snes_MixedFlowSolver converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
	ierr = DMDAGetCorners(ctx->daScal,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
	for (k = zs; k < zs+zm; k++) {
		for (j = ys; j < ys+ym; j++) {
			for (i = xs; i < xs+xm; i++) {
				Press_array[k][j][i] = VelnPress_array[k][j][i][3];
				for (c = 0; c < veldof; c++) {
					vel_array[k][j][i][c] =  VelnPress_array[k][j][i][c];
				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(ctx->daScal,fields->pressure,&Press_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,fields->velocity,&vel_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,fields->VelnPress,&VelnPress_array);CHKERRQ(ierr);
  ierr = VecMin(fields->velocity,PETSC_NULL,&Velmin);CHKERRQ(ierr);
	ierr = VecMax(fields->velocity,PETSC_NULL,&Velmax);CHKERRQ(ierr);
  ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
	ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Velocity min / max:     %e %e\n",Velmin,Velmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Pressure min / max:     %e %e\n",Pmin,Pmax);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESIJacobian"
extern PetscErrorCode FormSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
	PetscErrorCode    ierr;
	VFCtx             *ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	*str = DIFFERENT_NONZERO_PATTERN;
	ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KVelP,*Jac,*str);
	ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*Jac != *Jacpre) {
		ierr = MatAssemblyBegin(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}