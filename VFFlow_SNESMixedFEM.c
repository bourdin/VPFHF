/*
 VFFlow_MixedFEM.c
 A mixed finite elements Darcy solver based on the method in
 Masud, A. and Hughes, T. J. (2002). A stabilized mixed finite element method for
 Darcy flow. Computer Methods in Applied Mechanics and Engineering, 191(3940):43414370.
 
 (c) 2011-2012 C. Chukwudozie, LSU
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
/* #include "PetscFixes.h" */
#include "VFFlow_SNESMixedFEM.h"
#include "VFFlow_KSPMixedFEM.h"

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
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PetscOptionsBegin(PETSC_COMM_WORLD,PETSC_NULL,"","");CHKERRQ(ierr);
	{
		ctx->units    = UnitaryUnits;
		ierr          = PetscOptionsEnum("-flowunits","\n\tFlow solver","",FlowUnitName,(PetscEnum)ctx->units,(PetscEnum*)&ctx->units,PETSC_NULL);CHKERRQ(ierr);
			//	ctx->flowcase = ALLPRESSUREBC; 
		ctx->flowcase = ALLNORMALFLOWBC;
		ierr          = PetscOptionsEnum("-flow boundary conditions","\n\tFlow solver","",FlowBC_Case,(PetscEnum)ctx->flowcase,(PetscEnum*)&ctx->flowcase,PETSC_NULL);CHKERRQ(ierr);
	}
	ierr = PetscOptionsEnd();CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
	if (comm_size == 1) {
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->KVelPlhs);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATSEQAIJ,&ctx->JacVelP);CHKERRQ(ierr);
	} else {
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelP);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->KVelPlhs);CHKERRQ(ierr);
		ierr = DMCreateMatrix(ctx->daFlow,MATMPIAIJ,&ctx->JacVelP);CHKERRQ(ierr);
	}
	ierr = MatZeroEntries(ctx->JacVelP);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->KVelPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->JacVelP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->PreFlowFields);CHKERRQ(ierr);
	ierr = VecSet(ctx->PreFlowFields,0.0);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->FlowFunct);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->RHSVelP,"RHS vector of flow equation");CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->FlowFunct,"RHS of TS flow solver");CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSVelPpre,0.);CHKERRQ(ierr);
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesVelP);CHKERRQ(ierr);
	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
	ierr = SETFlowBC(&ctx->bcFlow[0],ctx->flowcase);CHKERRQ(ierr);
	ierr = SETSourceTerms(ctx->Source,ctx->flowprop);
	ierr = SETBoundaryTerms(ctx,fields);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MixedFEMSNESFlowSolverFinalize"
extern PetscErrorCode MixedFEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MatDestroy(&ctx->KVelP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KVelPlhs);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->JacVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PreFlowFields);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->FlowFunct);CHKERRQ(ierr);
	ierr = SNESDestroy(&ctx->snesVelP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSVelPpre);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MixedFEMSNESSolve"
extern PetscErrorCode MixedFlowFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	PetscViewer        viewer;
	KSPConvergedReason reason;	
	PetscReal          ****VelnPress_array;
	PetscReal          ***Press_array;
	PetscReal          ****vel_array;
	PetscInt           i,j,k,c,veldof = 3;
	PetscInt           xs,xm,ys,ym,zs,zm;
	PetscInt           its;
	
	PetscFunctionBegin;	
/*	temporary created permfield in ctx so permeability can be in ctx	*/
	ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Perm);CHKERRQ(ierr);
	ierr = VecSet(ctx->Perm,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->vfperm,ctx->Perm);CHKERRQ(ierr);
/*	temporary created bcfield in ctx so bc values can be in ctx	*/
	ierr = DMCreateGlobalVector(ctx->daFlow,&ctx->FlowBC);CHKERRQ(ierr);
	ierr = VecSet(ctx->FlowBC,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->FlowBCArray,ctx->FlowBC);CHKERRQ(ierr);
/* Copying previous time step solution to the previous time step vector	*/
	ierr = VecCopy(fields->VelnPress,ctx->PreFlowFields);CHKERRQ(ierr);
	ierr = SNESSetFunction(ctx->snesVelP,ctx->FlowFunct,FormSNESIFunction,ctx);CHKERRQ(ierr);
    ierr = SNESSetJacobian(ctx->snesVelP,ctx->JacVelP,ctx->JacVelP,FormSNESIJacobian,ctx);CHKERRQ(ierr);
	ierr = SNESMonitorSet(ctx->snesVelP,MixedFEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = SNESSolve(ctx->snesVelP,PETSC_NULL,fields->VelnPress);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,ctx->RHSVelPpre);CHKERRQ(ierr);
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
	ierr = VecDestroy(&ctx->FlowBC);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MixedFEMSNESMonitor"
extern PetscErrorCode MixedFEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;
	Vec				solution;
	
	PetscFunctionBegin;
	ierr = SNESGetSolution(snes,&solution);CHKERRQ(ierr);	
	ierr = VecNorm(solution,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(solution,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(solution,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"iter_step %D :solution norm = %G, max sol. value  = %G, min sol. value = %G\n",its,norm,vmax,vmin);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormSNESIJacobian"
extern PetscErrorCode FormSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
	PetscErrorCode ierr;
	VFCtx				*ctx=(VFCtx*)user;
	PetscViewer        viewer;
	
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

#undef __FUNCT__
#define __FUNCT__ "FormSNESMatricesnVector"
extern PetscErrorCode FormSNESMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l;
	PetscInt       veldof = 3;
	PetscInt       c;
	PetscReal      ****perm_array;
	PetscReal      ****coords_array;
	PetscReal      ****RHS_array;
	PetscReal      *RHS_local;
	Vec            RHS_localVec;
	Vec            perm_local;
	PetscReal      hx,hy,hz;
	PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local,*KS_local,*KDlhs_local;
	PetscReal      *KArhs_local,*KBrhs_local,*KDrhs_local,*KBTransrhs_local;
	PetscReal      beta_c,mu,gx,gy,gz;
	PetscReal	   theta,timestepsize;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row,*row1;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      M_inv;

	PetscFunctionBegin;
	M_inv     = ctx->flowprop.M_inv;
	beta_c = ctx->flowprop.beta;
	theta = ctx->flowprop.theta;
	timestepsize = ctx->flowprop.timestepsize;
	mu     = ctx->flowprop.mu;
	gx     = ctx->flowprop.g[0];
	gy     = ctx->flowprop.g[1];
	gz     = ctx->flowprop.g[2];
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
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KA_local,
						nrow*nrow,PetscReal,&KB_local,
						nrow*nrow,PetscReal,&KD_local,
						nrow*nrow,PetscReal,&KBTrans_local,
						nrow*nrow,PetscReal,&KS_local);CHKERRQ(ierr);
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KArhs_local,
						nrow*nrow,PetscReal,&KBrhs_local,
						nrow*nrow,PetscReal,&KBTransrhs_local,
						nrow*nrow,PetscReal,&KDrhs_local,
						nrow*nrow,PetscReal,&KDlhs_local);CHKERRQ(ierr);
	ierr = PetscMalloc3(nrow,PetscReal,&RHS_local,
						nrow,MatStencil,&row,
						nrow,MatStencil,&row1);CHKERRQ(ierr);	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
					//This computes the local contribution of the global A matrix
				ierr = FLow_MatA(KA_local,&ctx->e3D,ek,ej,ei);CHKERRQ(ierr);
				for (l = 0; l < nrow*nrow; l++) {
					KS_local[l] = 2*M_inv*KA_local[l];
					KArhs_local[l] = -1.*(1.-theta)*KA_local[l];
					KA_local[l] = theta*KA_local[l];
				}
				for (c = 0; c < veldof; c++) {
					ierr = FLow_MatB(KB_local,&ctx->e3D,ek,ej,ei,c);CHKERRQ(ierr);
					ierr = FLow_MatBTranspose(KBTrans_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array);CHKERRQ(ierr);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								row[l].i  = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = c;
								row1[l].i = ei+i;row1[l].j = ej+j;row1[l].k = ek+k;row1[l].c = 3;
							}
						}
					}
					for (l = 0; l < nrow*nrow; l++) {
						KBrhs_local[l] = -1.*timestepsize*(1.-theta)*KB_local[l];
						KBTransrhs_local[l] = -1.*(1.-theta)*KBTrans_local[l];
						KB_local[l] = timestepsize*theta*KB_local[l];
						KBTrans_local[l] = theta*KBTrans_local[l];
				}
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row1,nrow,row,KB_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(K,nrow,row,nrow,row1,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,KArhs_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(Klhs,nrow,row1,nrow,row,KBrhs_local,ADD_VALUES);CHKERRQ(ierr);
					ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row1,KBTransrhs_local,ADD_VALUES);CHKERRQ(ierr);
				}
				ierr = FLow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 3;
						}
					}
				}
				for (l = 0; l < nrow*nrow; l++) {
					KDlhs_local[l] = KS_local[l]+timestepsize*theta*KD_local[l];
					KDrhs_local[l] = KS_local[l]-timestepsize*(1.-theta)*KD_local[l];
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,KDlhs_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,KDrhs_local,ADD_VALUES);CHKERRQ(ierr);
					//Assembling the righthand side vector f
				for (c = 0; c < veldof; c++) {
					ierr = FLow_Vecf(RHS_local,&ctx->e3D,ek,ej,ei,c,ctx->flowprop,perm_array);CHKERRQ(ierr);
					for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
						for (j = 0; j < ctx->e3D.nphiy; j++) {
							for (i = 0; i < ctx->e3D.nphix; i++,l++) {
								RHS_array[ek+k][ej+j][ei+i][c] += RHS_local[l];
							}
						}
					}
				}
					//Assembling the righthand side vector g
				ierr = FLow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += timestepsize*RHS_local[l];
						}
					}
				}
				ierr = VecApplyWellFlowBC(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += timestepsize*RHS_local[l];
						}
					}
				}
			}
		}
	}	
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = ApplySNESJacobianBC(K,Klhs,&ctx->bcFlow[0]);CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = PetscFree5(KA_local,KB_local,KD_local,KBTrans_local,KS_local);CHKERRQ(ierr);
	ierr = PetscFree5(KArhs_local,KBrhs_local,KBTransrhs_local,KDrhs_local,KDlhs_local);CHKERRQ(ierr);
	ierr = PetscFree3(RHS_local,row,row1);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESIFunction"
extern PetscErrorCode FormSNESIFunction(SNES snes,Vec VelnPress,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	PetscViewer     viewer;
	Vec				vec;
	PetscReal			theta,one_minus_theta;

	
	PetscFunctionBegin;
	theta = ctx->flowprop.theta;
	one_minus_theta = (1.-theta);
	ierr = VecDuplicate(ctx->RHSVelP,&vec);CHKERRQ(ierr);
	ierr = FormSNESMatricesnVector(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,ctx);CHKERRQ(ierr);
	ierr = VecAXPBY(ctx->RHSVelP,one_minus_theta,theta,ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = VecApplyTSFlowBC(ctx->RHSVelP,ctx->FlowBC,&ctx->bcFlow[0],ctx);CHKERRQ(ierr);
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = VecSet(vec,0.0);CHKERRQ(ierr);
	ierr = MatMult(ctx->KVelP,VelnPress,Func);CHKERRQ(ierr);	
	ierr = MatMultAdd(ctx->KVelPlhs,ctx->PreFlowFields,ctx->RHSVelP,vec);CHKERRQ(ierr);	
	ierr = VecAXPY(Func,-1.0,vec);CHKERRQ(ierr);
	ierr = VecDestroy(&vec);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplySNESJacobianBC"
extern PetscErrorCode ApplySNESJacobianBC(Mat K,Mat Klhs,FLOWBC *BC)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       i,j,k,c;
	MatStencil    *row;
	PetscReal      one=1.;
	PetscInt       numBC=0,l=0;
	PetscInt       dim,dof;
	DM				da;
	PetscReal      zero=0.0;
	
	PetscFunctionBegin;
	ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   &dof,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	/*
	 Compute the number of boundary nodes on each processor. 
	 Edges and corners are counted multiple times (2 and 3 resp)
	 */
	for (c = 0; c < dof; c++){
		if (xs == 0       && BC[c].face[X0] != NOBC)             numBC += ym * zm;
		if (xs + xm == nx && BC[c].face[X1] != NOBC)             numBC += ym * zm;
		if (ys == 0       && BC[c].face[Y0] != NOBC)             numBC += xm * zm;
		if (ys + ym == ny && BC[c].face[Y1] != NOBC)             numBC += xm * zm;
		if (zs == 0       && BC[c].face[Z0] != NOBC && dim == 3) numBC += xm * ym;
		if (zs + zm == nz && BC[c].face[Z1] != NOBC && dim == 3) numBC += xm * ym;
		if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NOBC) numBC++;
		if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NOBC) numBC++;
		if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NOBC) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NOBC) numBC++;
		if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NOBC && dim == 3) numBC++;
		if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NOBC && dim == 3) numBC++;
		if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NOBC && dim == 3) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NOBC && dim == 3) numBC++;
	}
	ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
	/*
	 Create an array of rows to be zeroed out
	 */
	/*
	 i == 0
	 */
	for (c = 0; c < dof; c++) {
		if (xs == 0 && BC[c].face[X0] != NOBC) {
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
		if (xs + xm == nx && BC[c].face[X1] != NOBC) {
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
		if (ys == 0 && BC[c].face[Y0] != NOBC) {
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
		if (ys + ym == ny && BC[c].face[Y1] != NOBC) {
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
			if (zs == 0 && BC[c].face[Z0] != NOBC) {
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
			if (zs + zm == nz && BC[c].face[Z1] != NOBC) {
				for (j = ys; j < ys + ym; j++) {
					for (i = xs; i < xs + xm; i++) {
						row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
						l++;
					}
				}
			}
		}
		if (xs == 0       && ys == 0       && zs == 0       && BC[c].vertex[X0Y0Z0] != NOBC) { 
			row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys == 0       && zs + zm == nz && BC[c].vertex[X0Y0Z1] != NOBC && dim ==3) { 
			row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs == 0       && BC[c].vertex[X0Y1Z0] != NOBC) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs + zm == nz && BC[c].vertex[X0Y1Z1] != NOBC && dim ==3) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs == 0       && BC[c].vertex[X1Y0Z0] != NOBC) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs + zm == nz && BC[c].vertex[X1Y0Z1] != NOBC && dim ==3) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs == 0       && BC[c].vertex[X1Y1Z0] != NOBC) { 
			row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && BC[c].vertex[X1Y1Z1] != NOBC && dim ==3) { 
			row[l].i = nx=1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		
	}
	ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


