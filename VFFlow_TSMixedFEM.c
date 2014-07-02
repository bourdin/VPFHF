/*
 VFFlow_TSMixedFEM.c
 A mixed finite elements Darcy solver based on the method in
 Masud, A. and Hughes, T. J. (2002). A stabilized mixed finite element method for
 Darcy flow. Computer Methods in Applied Mechanics and Engineering, 191(3940):43414370.
 
 (c) 2011-2012 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_TSMixedFEM.h"
#include "VFFlow_KSPMixedFEM.h"

/*
 VFFlow_DarcyMixedFEMTS
*/
#undef __FUNCT__
#define __FUNCT__ "MixedFEMTSFlowSolverInitialize"
extern PetscErrorCode MixedFEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscErrorCode ierr;
  PetscFunctionBegin;

	ierr = TSCreate(PETSC_COMM_WORLD,&ctx->tsVelP);CHKERRQ(ierr);
//  ierr = TSAppendOptionsPrefix(ctx->tsVelP,"TSVelP_");CHKERRQ(ierr);
	ierr = TSSetDM(ctx->tsVelP,ctx->daFlow);CHKERRQ(ierr);
	ierr = TSSetProblemType(ctx->tsVelP,TS_LINEAR);CHKERRQ(ierr);
	ierr = TSSetType(ctx->tsVelP,TSBEULER);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMTSFlowSolverFinalize"
extern PetscErrorCode MixedFEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = TSDestroy(&ctx->tsVelP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFlowFEMTSSolve"
extern PetscErrorCode MixedFlowFEMTSSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	TSConvergedReason reason;	
	PetscReal          ****VelnPress_array;
	PetscReal          ***Press_array;
	PetscReal          ****vel_array;
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscReal           Velmin,Velmax;
	PetscReal           Pmin,Pmax;

	PetscFunctionBegin;
  /*	temporary created V in ctx so V-field in fields can be in ctx	*/
	ierr = VecSet(ctx->V,1.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->V,ctx->V);CHKERRQ(ierr);
  
  /*	temporary created permfield in ctx so permeability an be in ctx */
	ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Perm);CHKERRQ(ierr);
	ierr = VecSet(ctx->Perm,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->vfperm,ctx->Perm);CHKERRQ(ierr);
	ierr = TSSetIFunction(ctx->tsVelP,PETSC_NULL,FormIFunction,ctx);CHKERRQ(ierr);
  ierr = TSSetIJacobian(ctx->tsVelP,ctx->JacVelP,ctx->JacVelP,FormIJacobian,ctx);CHKERRQ(ierr);
	ierr = TSSetRHSFunction(ctx->tsVelP,PETSC_NULL,FormFunction,ctx);CHKERRQ(ierr);
	ierr = TSSetSolution(ctx->tsVelP,fields->VelnPress);CHKERRQ(ierr);
	ierr = TSSetInitialTimeStep(ctx->tsVelP,0.0,ctx->timevalue);CHKERRQ(ierr);
  ierr = TSSetDuration(ctx->tsVelP,ctx->maxtimestep,ctx->maxtimevalue);CHKERRQ(ierr);
//	ierr = FormInitialSolution(fields->VelnPress,fields->FlowBCArray,&ctx->bcFlow[0],ctx);CHKERRQ(ierr);
//	if (ctx->verbose > 1) {
		ierr = TSMonitorSet(ctx->tsVelP,MixedFEMTSMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
//	}	
	ierr = TSSetFromOptions(ctx->tsVelP);CHKERRQ(ierr);	
  ierr = TSSolve(ctx->tsVelP,fields->VelnPress);CHKERRQ(ierr);
  ierr = TSGetTimeStepNumber(ctx->tsVelP,&ctx->timestep);CHKERRQ(ierr);
	ierr = TSGetConvergedReason(ctx->tsVelP,&reason);CHKERRQ(ierr);
//	 PetscPrintf(PETSC_COMM_WORLD,"%s at time %G after %D steps\n",TSConvergedReasons[reason],ctx->timevalue,ctx->timestep);
/*
	if (reason < 0) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"[ERROR] tsVelP diverged with reason %d\n",(int)reason);CHKERRQ(ierr);
	} else {
		ierr = TSGetIterationNumber(ctx->snesVelP,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"      tsVelP converged in %d iterations %d.\n",(int)its,(int)reason);CHKERRQ(ierr);
	}
 */
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
	ierr = VecDestroy(&ctx->Perm);CHKERRQ(ierr);
  ierr = VecMin(fields->velocity,PETSC_NULL,&Velmin);CHKERRQ(ierr);
	ierr = VecMax(fields->velocity,PETSC_NULL,&Velmax);CHKERRQ(ierr);
  ierr = VecMin(fields->pressure,PETSC_NULL,&Pmin);CHKERRQ(ierr);
	ierr = VecMax(fields->pressure,PETSC_NULL,&Pmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Velocity min / max:     %e %e\n",Velmin,Velmax);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"      Pressure min / max:     %e %e\n",Pmin,Pmax);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MixedFEMTSMonitor"
extern PetscErrorCode MixedFEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec VelnPress,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;

	PetscFunctionBegin;
	ierr = VecNorm(VelnPress,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(VelnPress,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(VelnPress,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"timestep %D: time %G, solution norm %G, max %G, min %G\n",timestep,timevalue,norm,vmax,vmin);CHKERRQ(ierr);		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormIFunction"
extern PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = FormTSMatricesnVector(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,ctx);CHKERRQ(ierr);
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = MatMult(ctx->KVelP,VelnPress,Func);CHKERRQ(ierr);	
	ierr = MatMultAdd(ctx->KVelPlhs,VelnPressdot,Func,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunction"
extern PetscErrorCode FormFunction(TS ts,PetscReal t,Vec vec1,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormIJacobian"
extern PetscErrorCode FormIJacobian(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,PetscReal shift,Mat Jac,Mat Jacpre,void *user)
{
	PetscErrorCode ierr;
	VFCtx				*ctx=(VFCtx*)user;

	PetscFunctionBegin;
	ierr = MatCopy(ctx->KVelP,Jacpre,DIFFERENT_NONZERO_PATTERN);
	ierr = MatAXPY(Jacpre,shift,ctx->KVelPlhs,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (Jac != Jacpre) {
    ierr = MatCopy(Jacpre,Jac,DIFFERENT_NONZERO_PATTERN);
		ierr = MatAssemblyBegin(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "MatApplyTSVelocityBC"
extern PetscErrorCode MatApplyTSVelocityBC(Mat K,Mat Klhs,VFBC *bcQ)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       i,j,k,c;
	MatStencil    *row;
	PetscReal      one=1.;
	PetscInt       numBC=0,l=0;
	PetscInt       dim,dof=3;
	DM				da;
	PetscReal      zero=0.0;
	
	PetscFunctionBegin;
  ierr = MatGetDM(K,&da);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	for (c = 0; c < dof; c++){
		if (xs == 0       && bcQ[c].face[X0] == FIXED)             numBC += ym * zm;
		if (xs + xm == nx && bcQ[c].face[X1] == FIXED)             numBC += ym * zm;
		if (ys == 0       && bcQ[c].face[Y0] == FIXED)             numBC += xm * zm;
		if (ys + ym == ny && bcQ[c].face[Y1] == FIXED)             numBC += xm * zm;
		if (zs == 0       && bcQ[c].face[Z0] == FIXED && dim == 3) numBC += xm * ym;
		if (zs + zm == nz && bcQ[c].face[Z1] == FIXED && dim == 3) numBC += xm * ym;
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == FIXED) numBC++;
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == FIXED) numBC++;
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == FIXED) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == FIXED) numBC++;
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == FIXED && dim == 3) numBC++;
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == FIXED && dim == 3) numBC++;
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == FIXED && dim == 3) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == FIXED && dim == 3) numBC++;
	}
	ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
	/*
	 Create an array of rows to be zeroed out
	 */
	/*
	 i == 0
	 */
	for (c = 0; c < dof; c++) {
		if (xs == 0 && bcQ[c].face[X0] == FIXED) {
			for (k = zs; k < zs + zm; k++) {
				for (j = ys; j < ys + ym; j++) {
					row[l].i = 0; row[l].j = j; row[l].k = k; row[l].c = c; 
					l++;
				}
			}
		}
		/* 
		 i == nx-1
		 */
		if (xs + xm == nx && bcQ[c].face[X1] == FIXED) {
			for (k = zs; k < zs + zm; k++) {
				for (j = ys; j < ys + ym; j++) {
					row[l].i = nx-1; row[l].j = j; row[l].k = k; row[l].c = c; 
					l++;
				}
			}
		}
		/*
		 y == 0
		 */
		if (ys == 0 && bcQ[c].face[Y0] == FIXED) {
			for (k = zs; k < zs + zm; k++) {
				for (i = xs; i < xs + xm; i++) {
					row[l].i = i; row[l].j = 0; row[l].k = k; row[l].c = c; 
					l++;
				}
			}
		}
		/*
		 y == ny-1
		 */
		if (ys + ym == ny && bcQ[c].face[Y1] == FIXED) {
			for (k = zs; k < zs + zm; k++) {
				for (i = xs; i < xs + xm; i++) {
					row[l].i = i; row[l].j = ny-1; row[l].k = k; row[l].c = c; 
					l++;
				}
			}
		}
		if (dim==3){
			/*
			 z == 0
			 */
			if (zs == 0 && bcQ[c].face[Z0] == FIXED) {
				for (j = ys; j < ys + ym; j++) {
					for (i = xs; i < xs + xm; i++) {
						row[l].i = i; row[l].j = j; row[l].k = 0; row[l].c = c; 
						l++;
					}
				}
			}
			/*
			 z == nz-1
			 */
			if (zs + zm == nz && bcQ[c].face[Z1] == FIXED) {
				for (j = ys; j < ys + ym; j++) {
					for (i = xs; i < xs + xm; i++) {
						row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
						l++;
					}
				}
			}
		}
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == FIXED) { 
			row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == FIXED && dim ==3) { 
			row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == FIXED) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == FIXED && dim ==3) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == FIXED) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == FIXED && dim ==3) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == FIXED) { 
			row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == FIXED && dim ==3) { 
			row[l].i = nx=1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		
	}
	ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormTSMatricesnVector"
extern PetscErrorCode FormTSMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l,ii;
	PetscInt       veldof = 3;
	PetscInt       c;
	PetscReal      ****perm_array;
	PetscReal      ****coords_array;
	PetscReal      ****RHS_array;
	PetscReal      *RHS_local;
	Vec            RHS_localVec;
	Vec            perm_local;
	PetscReal      hx,hy,hz;
	PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local,*Klhs_local;
	PetscReal      alpha_c;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row,*row1;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      M_inv;
	FACE           face;
	PetscReal      ***prebc_array;
	Vec            prebc_local;
	PetscReal      hwx, hwy, hwz;	
	PetscReal      ***v_array;
	Vec            v_local;

	PetscFunctionBegin;
	M_inv     = ctx->flowprop.M_inv;
	alpha_c = ctx->flowprop.alpha;
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
		// This line ensures that the number of cells is one less than the number of nodes. Force processing of cells to stop once the second to the last node is processed 
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);
	ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&prebc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->PresBCArray,INSERT_VALUES,prebc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->PresBCArray,INSERT_VALUES,prebc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,prebc_local,&prebc_array);CHKERRQ(ierr); 
  ierr = DMGetLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->V,INSERT_VALUES,v_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = PetscMalloc5(nrow*nrow,&KA_local,
						nrow*nrow,&KB_local,
						nrow*nrow,&KD_local,
						nrow*nrow,&KBTrans_local,
						nrow*nrow,&Klhs_local);CHKERRQ(ierr);
	ierr = PetscMalloc3(nrow,&RHS_local,
						nrow,&row,
						nrow,&row1);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				/*	This computes the local contribution of the global A matrix	*/
        ierr = VF_MatA_local(Klhs_local,&ctx->e3D,ek,ej,ei,v_array);CHKERRQ(ierr);
				for (l = 0; l < nrow*nrow; l++) {
					Klhs_local[l] = -1.*M_inv*Klhs_local[l]/alpha_c;
				}
				for (c = 0; c < veldof; c++) {
          ierr = Flow_MatA(KA_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
					ierr = Flow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c,v_array);CHKERRQ(ierr);
					ierr = Flow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,v_array);CHKERRQ(ierr);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								row[l].i  = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = c;
								row1[l].i = ei+i;row1[l].j = ej+j;row1[l].k = ek+k;row1[l].c = 3;
							}
						}
					}
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KB_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
				}
				ierr = Flow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(K,nrow,row1,nrow,row1,KD_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row1,nrow,row1,Klhs_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the righthand side vector f*/
				for (c = 0; c < veldof; c++) {
					ierr = Flow_Vecf(RHS_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,v_array);CHKERRQ(ierr);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
							}
						}
					}
				}
				/*Assembling the right hand side vector g*/
				ierr = Flow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
						}
					}
				}
				if(ctx->hasFluidSources){
					ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx,v_array);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
							}
						}
					}
				}
				if (ei == 0) {
					/*					 Face X0			*/
					face = X0;	
					ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i][0] -= RHS_local[l];
								}
							}
						}
					}
				}
				if (ei == nx-1) {
					/*					 Face X1		*/
					face = X1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphix; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphiz; i++, l++) {
									RHS_array[ek+k][ej+j][ei+1][0] += RHS_local[l];
								}
							}
						}
					}
				}				
				if (ej == 0) {
					/*					 Face Y0		*/
					face = Y0;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i][1] -= RHS_local[l];
								}
							}
						}
					}
				}
				if (ej == ny-1) {
					/*					 Face Y1		*/
					face = Y1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiy; k++){
							for (j = 0; j < ctx->e2D.nphiz; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+1][ei+i][1] += RHS_local[l];
								}
							}
						}
					}
				}
				if (ek == 0) {
					/*					 Face Z0		*/
					face = Z0;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+k][ej+j][ei+i][2] -= RHS_local[l];
								}
							}
						}
					}
				}
				if (ek == nz-1) {
					/*					 Face Z1		*/
					face = Z1;
					ierr = VFCartFEElement2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == FIXED) {
						ierr = VecApplyPressureBC(RHS_local,prebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array,v_array);CHKERRQ(ierr);
						for (l=0,k = 0; k < ctx->e2D.nphiz; k++){
							for (j = 0; j < ctx->e2D.nphiy; j++) {
								for (i = 0; i < ctx->e2D.nphix; i++, l++) {
									RHS_array[ek+1][ej+j][ei+i][2] += RHS_local[l];
								}
							}
						}
					}
				}
			}
		}
	}		
	if(ctx->hasFlowWells){
		for(ii = 0; ii < ctx->numWells; ii++){
			for (ek = zs; ek < zs+zm; ek++) {
				for (ej = ys; ej < ys+ym; ej++) {
					for (ei = xs; ei < xs+xm; ei++) {
						if(  
						   ((coords_array[ek][ej][ei+1][0] >= ctx->well[ii].coords[0]) && (coords_array[ek][ej][ei][0] <= ctx->well[ii].coords[0] ))	&&
						   ((coords_array[ek][ej+1][ei][1] >= ctx->well[ii].coords[1]) && (coords_array[ek][ej][ei][1] <= ctx->well[ii].coords[1] ))	&&
						   ((coords_array[ek+1][ej][ei][2] >= ctx->well[ii].coords[2]) && (coords_array[ek][ej][ei][2] <= ctx->well[ii].coords[2] ))
						   )
						{
              hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
              hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
              hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
              ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
							hwx = (ctx->well[ii].coords[0]-coords_array[ek][ej][ei][0])/hx;
							hwy = (ctx->well[ii].coords[1]-coords_array[ek][ej][ei][1])/hy;
							hwz = (ctx->well[ii].coords[2]-coords_array[ek][ej][ei][2])/hz;
							if(ctx->well[ii].condition == RATE){
								ierr = VecApplyWellFlowRate(RHS_local,&ctx->e3D,ctx->well[ii].Qw,hwx,hwy,hwz,ek,ej,ei,v_array);CHKERRQ(ierr);
								for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
									for (j = 0; j < ctx->e3D.nphiy; j++) {
										for (i = 0; i < ctx->e3D.nphix; i++,l++) {
											if(ctx->well[ii].type == INJECTOR){
												RHS_array[ek+k][ej+j][ei+i][3] -= RHS_local[l];
											}
											else if(ctx->well[ii].type == PRODUCER)
											{
												RHS_array[ek+k][ej+j][ei+i][3] -= -1*RHS_local[l];
											}
										}
									}
								}
							}
							else if(ctx->well[ii].condition == PRESSURE){
								
							}
							goto outer;
						}
					}
				}
			}
		outer:;
		}
	}
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatApplyTSVelocityBC(K,Klhs,&ctx->bcQ[0]);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daFlow,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daFlow,&RHS_localVec);CHKERRQ(ierr);	
	ierr = VecApplyTSVelocityBC(RHS,ctx->VelBCArray,&ctx->bcQ[0],ctx);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,prebc_local,&prebc_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daScal,&prebc_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,v_local,&v_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&v_local);CHKERRQ(ierr);
	ierr = PetscFree5(KA_local,KB_local,KD_local,KBTrans_local,Klhs_local);CHKERRQ(ierr);
	ierr = PetscFree3(RHS_local,row,row1);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplyTSVelocityBC"
extern PetscErrorCode VecApplyTSVelocityBC(Vec RHS,Vec BCV, VFBC *bcQ,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof=3;
	PetscInt       i,j,k,c;
	PetscReal		****velbc_array;
	PetscReal		****RHS_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daFlow,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,BCV,&velbc_array);CHKERRQ(ierr); 	
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr); 	

	for (c = 0; c < dof; c++) {
		if (xs == 0) {
			i = 0;
			if (bcQ[c].face[X0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs+xm == nx) {
			i = nx-1;
			if (bcQ[c].face[X1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys == 0) {
			j = 0;
			if (bcQ[c].face[Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (ys+ym == ny) {
			j = ny-1;
			if (bcQ[c].face[Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
						
					}
				}
			}
		}
		if (zs == 0) {
			k = 0;
			if (bcQ[c].face[Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (zs+zm == nz) {
			k = nz-1;
			if (bcQ[c].face[Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
					}
				}
			}
		}
		if (xs == 0 && zs == 0) {
			k = 0;i = 0;
			if (bcQ[c].edge[X0Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs == 0) {
			k = 0;i = nx-1;
			if (bcQ[c].edge[X1Z0] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y0Z0] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y1Z0] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1;i = 0;
			if (bcQ[c].edge[X0Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1;i = nx-1;
			if (bcQ[c].edge[X1Z1] == FIXED) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;
			if (bcQ[c].edge[Y0Z1] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;
			if (bcQ[c].edge[Y1Z1] == FIXED) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0) {
			j = 0;i = 0;
			if (bcQ[c].edge[X0Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys+ym == ny) {
			j = ny-1;i = 0;
			if (bcQ[c].edge[X0Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys == 0) {
			j = 0;i = nx-1;
			if (bcQ[c].edge[X1Y0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1;i = nx-1;
			if (bcQ[c].edge[X1Y1] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z0] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z1] == FIXED) {
				RHS_array[k][j][i][c] = velbc_array[k][j][i][c];
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,BCV,&velbc_array);CHKERRQ(ierr); 	
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr); 	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormInitialSolution"
extern PetscErrorCode FormInitialSolution(Vec VelnPress,Vec VelnPressBV, VFBC *bcP,VFBC *bcQ, VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof=3;
	PetscInt       i,j,k,c;
	DM             da;
	PetscReal		****UnPre_array;
	PetscReal      ****IniSol_array;
	
	PetscFunctionBegin;
  ierr = VecGetDM(VelnPress,&da);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,VelnPress,&IniSol_array);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,VelnPressBV,&UnPre_array);CHKERRQ(ierr);
		if (xs == 0) {
			i = 0;
			if (bcP[0].face[X0] == FIXED) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
					}
				}
			}
			else {
				for(c = 0; c < dof; c++){
					if (bcQ[c].face[X0] == FIXED) {
						for (k = zs; k < zs+zm; k++) {
							for (j = ys; j < ys+ym; j++) {
								IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
							}
						}
					}
				
				}
			
			}
		}

	
	if (xs+xm == nx) {
		i = nx-1;
		if (bcP[0].face[X1] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				for (j = ys; j < ys+ym; j++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].face[X1] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						for (j = ys; j < ys+ym; j++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}
	}

	if (ys == 0) {
		j = 0;
		if (bcP[0].face[Y0] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				for (i = xs; i < xs+xm; i++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].face[Y0] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						for (i = xs; i < xs+xm; i++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}
	}
	if (ys+ym == ny) {
		j = ny-1;
		if (bcP[0].face[Y1] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				for (i = xs; i < xs+xm; i++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].face[Y1] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						for (i = xs; i < xs+xm; i++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}
	}
	
	if (zs == 0) {
		k = 0;
		if (bcP[0].face[Z0] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].face[Z0] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}
	}
	if (zs+zm == nz) {
		k = nz-1;
		if (bcP[0].face[Z1] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
				for (i = xs; i < xs+xm; i++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].face[Z1] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
						for (i = xs; i < xs+xm; i++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}
	}
	
	if (xs == 0 && zs == 0) {
		k = 0;i = 0;
		if (bcP[0].edge[X0Z0] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
					IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
				}
			}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X0Z0] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
							IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
						}
					}
				}
				
			}
			
		}

	if (xs+xm == nx && zs == 0) {
		k = 0;i = nx-1;
		if (bcP[0].edge[X1Z0] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X1Z0] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}

		//
	if (ys == 0 && zs == 0) {
		k = 0;j = 0;
		if (bcP[0].edge[Y0Z0] == FIXED) {
			for (i = xs; i < xs+xm; i++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[Y0Z0] == FIXED) {
					for (i = xs; i < xs+xm; i++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
	
	if (ys+ym == ny && zs == 0) {
		k = 0;j = ny-1;
		if (bcP[0].edge[Y1Z0] == FIXED) {
			for (i = xs; i < xs+xm; i++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[Y1Z0] == FIXED) {
					for (i = xs; i < xs+xm; i++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
		//
	
	if (xs == 0 && zs+zm == nz) {
		k = nz-1;i = 0;
		if (bcP[0].edge[X0Z1] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X0Z1] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
	
	if (xs+xm == nx && zs+zm == nz) {
		k = nz-1;i = nx-1;
		if (bcP[0].edge[X1Z1] == FIXED) {
			for (j = ys; j < ys+ym; j++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X1Z1] == FIXED) {
					for (j = ys; j < ys+ym; j++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
		}
	}
		//
	if (ys == 0 && zs+zm == nz) {
		k = nz-1;j = 0;
		if (bcP[0].edge[Y0Z1] == FIXED) {
			for (i = xs; i < xs+xm; i++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[Y0Z1] == FIXED) {
					for (i = xs; i < xs+xm; i++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
	
	if (ys+ym == ny && zs+zm == nz) {
		k = nz-1;j = ny-1;
		if (bcP[0].edge[Y1Z1] == FIXED) {
			for (i = xs; i < xs+xm; i++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[Y1Z1] == FIXED) {
					for (i = xs; i < xs+xm; i++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}

		//
	if (xs == 0 && ys == 0) {
		j = 0;i = 0;
		if (bcP[0].edge[X0Y0] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X0Y0] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
	
	if (xs == 0 && ys+ym == ny) {
		j = ny-1;i = 0;
		if (bcP[0].edge[X0Y1] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X0Y1] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}

		//
	if (xs+xm == nx && ys == 0) {
		j = 0;i = nx-1;
		if (bcP[0].edge[X1Y0] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X1Y0] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}
	
	if (xs+xm == nx && ys+ym == ny) {
		j = ny-1;i = nx-1;
		if (bcP[0].edge[X1Y1] == FIXED) {
			for (k = zs; k < zs+zm; k++) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		}
		
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].edge[X1Y1] == FIXED) {
					for (k = zs; k < zs+zm; k++) {
						IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
			
		}
		
	}

		//
	
	if (xs == 0 && ys == 0 && zs == 0) {
		k = 0;j = 0;i = 0;
		if (bcP[0].vertex[X0Y0Z0] == FIXED) {
				IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
			}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X0Y0Z0] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
					}
				}
			}
	}
	if (xs+xm == nx && ys == 0 && zs == 0) {
		k = 0;j = 0;i = nx-1;
		if (bcP[0].vertex[X1Y0Z0] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X1Y0Z0] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}
		//
	if (xs == 0 && ys+ym == ny && zs == 0) {
		k = 0;j = ny-1;i = 0;
		if (bcP[0].vertex[X0Y1Z0] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X0Y1Z0] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}
	if (xs+xm == nx && ys+ym == ny && zs == 0) {
		k = 0;j = ny-1;i = nx-1;
		if (bcP[0].vertex[X1Y1Z0] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X1Y1Z0] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}

		//
	if (xs == 0 && ys == 0 && zs+zm == nz) {
		k = nz-1;j = 0;i = 0;
		if (bcP[0].vertex[X0Y0Z1] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X0Y0Z1] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}
	if (xs+xm == nx && ys == 0 && zs+zm == nz) {
		k = nz-1;j = 0;i = nx-1;
		if (bcP[0].vertex[X1Y0Z1] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X1Y0Z1] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}

		//
	
	if (xs == 0 && ys+ym == ny && zs+zm == nz) {
		k = nz-1;j = ny-1;i = 0;
		if (bcP[0].vertex[X0Y1Z1] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X0Y1Z1] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}
	if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
		k = nz-1;j = ny-1;i = nx-1;
		if (bcP[0].vertex[X1Y1Z1] == FIXED) {
			IniSol_array[k][j][i][3] = UnPre_array[k][j][i][3];
		}
		else {
			for(c = 0; c < dof; c++){
				if (bcQ[c].vertex[X1Y1Z1] == FIXED) {
					IniSol_array[k][j][i][c] = UnPre_array[k][j][i][c];
				}
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(da,VelnPress,&IniSol_array);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(da,VelnPressBV,&UnPre_array);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}







