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
	
	ierr = BCPInit(&ctx->bcP[0],ctx);
	ierr = BCQInit(&ctx->bcQ[0],ctx);
	
	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
		//	ierr = SETFlowBC(&ctx->bcFlow[0],ctx->flowcase);CHKERRQ(ierr);
	ierr = SETFlowBC(&ctx->bcP[0],&ctx->bcQ[0],ctx->flowcase);CHKERRQ(ierr);	
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
	PetscReal      *KA_local,*KB_local,*KD_local,*KBTrans_local,*KS_local;
	PetscReal      beta_c,alpha_c,mu,gx,gy,gz;
	PetscReal	   theta,timestepsize;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row,*row1;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      M_inv;
	Mat			   K1,K2;
	PetscReal		time_theta,time_one_minus_theta;
	FACE			face;
	PetscReal		****velnprebc_array;
	Vec				velnprebc_local;
	MatStructure	flg;
	PetscReal		hwx, hwy, hwz;	
	
	PetscFunctionBegin;
	flg = SAME_NONZERO_PATTERN;
	M_inv     = ctx->flowprop.M_inv;
	beta_c = ctx->flowprop.beta;
	alpha_c = ctx->flowprop.alpha;
	theta = ctx->flowprop.theta;
	timestepsize = ctx->flowprop.timestepsize;
	time_theta = theta * timestepsize;
	time_one_minus_theta = -1.*(1.-theta)*timestepsize;
	mu     = ctx->flowprop.mu;
	gx     = ctx->flowprop.g[0];
	gy     = ctx->flowprop.g[1];
	gz     = ctx->flowprop.g[2];
	
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
		// This line ensures that the number of cells is one less than the number of nodes. Force processing of cells to stop once the second to the last node is processed 
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
	ierr = MatDuplicate(K,MAT_COPY_VALUES,&K1);CHKERRQ(ierr);
	ierr = MatDuplicate(K,MAT_COPY_VALUES,&K2);CHKERRQ(ierr);
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
	
	ierr = DMGetLocalVector(ctx->daFlow,&velnprebc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daFlow,ctx->FlowBC,INSERT_VALUES,velnprebc_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daFlow,ctx->FlowBC,INSERT_VALUES,velnprebc_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,velnprebc_local,&velnprebc_array);CHKERRQ(ierr); 
	
	ierr = PetscMalloc5(nrow*nrow,PetscReal,&KA_local,
						nrow*nrow,PetscReal,&KB_local,
						nrow*nrow,PetscReal,&KD_local,
						nrow*nrow,PetscReal,&KBTrans_local,
						nrow*nrow,PetscReal,&KS_local);CHKERRQ(ierr);
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
					KS_local[l] = -2*M_inv*KA_local[l]/alpha_c;
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
					ierr = MatSetValuesStencil(K1,nrow,row,nrow,row,KA_local,ADD_VALUES);CHKERRQ(ierr);			
					ierr = MatSetValuesStencil(K1,nrow,row1,nrow,row,KB_local,ADD_VALUES);CHKERRQ(ierr);		
					ierr = MatSetValuesStencil(K1,nrow,row,nrow,row1,KBTrans_local,ADD_VALUES);CHKERRQ(ierr);
				}
				ierr = FLow_MatD(KD_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 3;
						}
					}
				}
				ierr = MatSetValuesStencil(K1,nrow,row,nrow,row,KD_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(K2,nrow,row,nrow,row,KS_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the righthand side vector f*/
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
				/*Assembling the righthand side vector g*/
				ierr = FLow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
						}
					}
				}
				if(ctx->hasFluidSources){
					ierr = VecApplyWellFlowBC(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hz,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hz);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
					ierr = CartFE_Element2DInit(&ctx->e2D,hx,hy);CHKERRQ(ierr);
					if (ctx->bcP[0].face[face] == VALUE) {
						ierr = VecApplyPressureBC(RHS_local,velnprebc_array,ek,ej,ei,face,&ctx->e2D,ctx->flowprop,perm_array);CHKERRQ(ierr);
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
								hwx = (ctx->well[ii].coords[0]-coords_array[ek][ej][ei][0])/(coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0]); 
								hwy = (ctx->well[ii].coords[1]-coords_array[ek][ej][ei][1])/(coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1]); 
								hwz = (ctx->well[ii].coords[2]-coords_array[ek][ej][ei][2])/(coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2]); 
								if(ctx->well[ii].condition == RATE){
									ierr = VecApplyWEllFlowRate(RHS_local,ctx->well[ii].Qw,hwx,hwy,hwz);CHKERRQ(ierr);
									for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
										for (j = 0; j < ctx->e3D.nphiy; j++) {
											for (i = 0; i < ctx->e3D.nphix; i++,l++) {
												if(ctx->well[ii].type == INJECTOR){
													RHS_array[ek+k][ej+j][ei+i][3] += -1*RHS_local[l];
												}
												else if(ctx->well[ii].type == PRODUCER)
												{
													RHS_array[ek+k][ej+j][ei+i][3] += RHS_local[l];
												}
											}
										}
									}
								}
								else if(ctx->well[ii].type == PRESSURE){
									
								}
								goto outer;
							}
					}
				}
			}
			outer:;
		}
	}
	ierr = MatAssemblyBegin(K2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAXPY(K,time_theta,K1,flg);CHKERRQ(ierr);
	ierr = MatAXPY(K,1.0,K2,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Klhs,time_one_minus_theta,K1,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Klhs,1.0,K2,flg);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatApplySNESVelocityBC(K,Klhs,&ctx->bcQ[0]);CHKERRQ(ierr);
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
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,velnprebc_local,&velnprebc_array);CHKERRQ(ierr); 
	ierr = DMRestoreLocalVector(ctx->daFlow,&velnprebc_local);CHKERRQ(ierr);
	ierr = PetscFree5(KA_local,KB_local,KD_local,KBTrans_local,KS_local);CHKERRQ(ierr);
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
	Vec             VecRHS;
	PetscReal		theta,timestepsize;
	PetscReal		dt_dot_theta,dt_dot_one_minus_theta;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	timestepsize = ctx->flowprop.timestepsize;
	theta = ctx->flowprop.theta;
	dt_dot_theta = timestepsize*theta;
	dt_dot_one_minus_theta = timestepsize*(1.-theta);
	ierr = VecDuplicate(ctx->RHSVelP,&VecRHS);CHKERRQ(ierr);
	ierr = FormSNESMatricesnVector(ctx->KVelP,ctx->KVelPlhs,ctx->RHSVelP,ctx);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSVelP,VecRHS);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,dt_dot_one_minus_theta,dt_dot_theta,ctx->RHSVelPpre);CHKERRQ(ierr);
	ierr = MatMultAdd(ctx->KVelPlhs,ctx->PreFlowFields,VecRHS,VecRHS);CHKERRQ(ierr);	
	ierr = VecApplySNESVelocityBC(VecRHS,ctx->FlowBC,&ctx->bcQ[0],ctx);CHKERRQ(ierr);
	ierr = MatMult(ctx->KVelP,VelnPress,Func);CHKERRQ(ierr);	
	ierr = VecAXPY(Func,-1.0,VecRHS);CHKERRQ(ierr);
	ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecApplySNESVelocityBC"
extern PetscErrorCode VecApplySNESVelocityBC(Vec RHS,Vec BCV, BC *bcQ,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       dim,dof=3;
	PetscInt       i,j,k,c;
	PetscReal		****bcv_array;
	PetscReal		****RHS_array;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daFlow,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daFlow,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,BCV,&bcv_array);CHKERRQ(ierr); 	
	ierr = DMDAVecGetArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr); 		
	for (c = 0; c < dof; c++) {
		if (xs == 0) {
			i = 0;
			if (bcQ[c].face[X0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
					}
				}
			}
		}
		if (xs+xm == nx) {
			i = nx-1;
			if (bcQ[c].face[X1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (j = ys; j < ys+ym; j++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
					}
				}
			}
		}
		if (ys == 0) {
			j = 0;
			if (bcQ[c].face[Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
					}
				}
			}
		}
		if (ys+ym == ny) {
			j = ny-1;
			if (bcQ[c].face[Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
						
					}
				}
			}
		}
		if (zs == 0) {
			k = 0;
			if (bcQ[c].face[Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
					}
				}
			}
		}
		if (zs+zm == nz) {
			k = nz-1;
			if (bcQ[c].face[Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					for (i = xs; i < xs+xm; i++) {
						RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
					}
				}
			}
		}
		if (xs == 0 && zs == 0) {
			k = 0;i = 0;
			if (bcQ[c].edge[X0Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs == 0) {
			k = 0;i = nx-1;
			if (bcQ[c].edge[X1Z0] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y0Z0] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs == 0) {
			k = 0;j = 0;
			if (bcQ[c].edge[Y1Z0] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && zs+zm == nz) {
			k = nz-1;i = 0;
			if (bcQ[c].edge[X0Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && zs+zm == nz) {
			k = nz-1;i = nx-1;
			if (bcQ[c].edge[X1Z1] == VALUE) {
				for (j = ys; j < ys+ym; j++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;
			if (bcQ[c].edge[Y0Z1] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;
			if (bcQ[c].edge[Y1Z1] == VALUE) {
				for (i = xs; i < xs+xm; i++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0) {
			j = 0;i = 0;
			if (bcQ[c].edge[X0Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys+ym == ny) {
			j = ny-1;i = 0;
			if (bcQ[c].edge[X0Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys == 0) {
			j = 0;i = nx-1;
			if (bcQ[c].edge[X1Y0] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs+xm == nx && ys+ym == ny) {
			j = ny-1;i = nx-1;
			if (bcQ[c].edge[X1Y1] == VALUE) {
				for (k = zs; k < zs+zm; k++) {
					RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
				}
			}
		}
		if (xs == 0 && ys == 0 && zs == 0) {
			k = 0;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z0] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs == 0) {
			k = 0;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z0] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = 0;
			if (bcQ[c].vertex[X0Y1Z0] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs == 0) {
			k = 0;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z0] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = 0;
			if (bcQ[c].vertex[X0Y0Z1] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys == 0 && zs+zm == nz) {
			k = nz-1;j = 0;i = nx-1;
			if (bcQ[c].vertex[X1Y0Z1] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs == 0 && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = 0;
			
			if (bcQ[c].vertex[X0Y1Z1] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
		if (xs+xm == nx && ys+ym == ny && zs+zm == nz) {
			k = nz-1;j = ny-1;i = nx-1;
			if (bcQ[c].vertex[X1Y1Z1] == VALUE) {
				RHS_array[k][j][i][c] = bcv_array[k][j][i][c];
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,BCV,&bcv_array);CHKERRQ(ierr); 	
	ierr = DMDAVecRestoreArrayDOF(ctx->daFlow,RHS,&RHS_array);CHKERRQ(ierr); 	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatApplySNESVelocityBC"
extern PetscErrorCode MatApplySNESVelocityBC(Mat K,Mat Klhs,BC *bcQ)
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
	ierr = PetscObjectQuery((PetscObject) K,"DM",(PetscObject *) &da); CHKERRQ(ierr);
	if (!da) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG," Matrix not generated from a DA");
	
	ierr = DMDAGetInfo(da,&dim,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,
					   PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	/*
	 Compute the number of boundary nodes on each processor. 
	 Edges and corners are counted multiple times (2 and 3 resp)
	 */
	for (c = 0; c < dof; c++){
		if (xs == 0       && bcQ[c].face[X0] == VALUE)             numBC += ym * zm;
		if (xs + xm == nx && bcQ[c].face[X1] == VALUE)             numBC += ym * zm;
		if (ys == 0       && bcQ[c].face[Y0] == VALUE)             numBC += xm * zm;
		if (ys + ym == ny && bcQ[c].face[Y1] == VALUE)             numBC += xm * zm;
		if (zs == 0       && bcQ[c].face[Z0] == VALUE && dim == 3) numBC += xm * ym;
		if (zs + zm == nz && bcQ[c].face[Z1] == VALUE && dim == 3) numBC += xm * ym;
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == VALUE) numBC++;
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == VALUE) numBC++;
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == VALUE) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == VALUE) numBC++;
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == VALUE && dim == 3) numBC++;
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == VALUE && dim == 3) numBC++;
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == VALUE && dim == 3) numBC++;
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == VALUE && dim == 3) numBC++;
	}
	ierr = PetscMalloc(numBC * sizeof(MatStencil),&row);CHKERRQ(ierr);
	/*
	 Create an array of rows to be zeroed out
	 */
	/*
	 i == 0
	 */
	for (c = 0; c < dof; c++) {
		if (xs == 0 && bcQ[c].face[X0] == VALUE) {
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
		if (xs + xm == nx && bcQ[c].face[X1] == VALUE) {
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
		if (ys == 0 && bcQ[c].face[Y0] == VALUE) {
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
		if (ys + ym == ny && bcQ[c].face[Y1] == VALUE) {
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
			if (zs == 0 && bcQ[c].face[Z0] == VALUE) {
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
			if (zs + zm == nz && bcQ[c].face[Z1] == VALUE) {
				for (j = ys; j < ys + ym; j++) {
					for (i = xs; i < xs + xm; i++) {
						row[l].i = i; row[l].j = j; row[l].k = nz-1; row[l].c = c; 
						l++;
					}
				}
			}
		}
		if (xs == 0       && ys == 0       && zs == 0       && bcQ[c].vertex[X0Y0Z0] == VALUE) { 
			row[l].i = 0; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys == 0       && zs + zm == nz && bcQ[c].vertex[X0Y0Z1] == VALUE && dim ==3) { 
			row[l].i = 0; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs == 0       && bcQ[c].vertex[X0Y1Z0] == VALUE) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs == 0       && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X0Y1Z1] == VALUE && dim ==3) { 
			row[l].i = 0; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs == 0       && bcQ[c].vertex[X1Y0Z0] == VALUE) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys == 0       && zs + zm == nz && bcQ[c].vertex[X1Y0Z1] == VALUE && dim ==3) { 
			row[l].i = nx-1; row[l].j = 0; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs == 0       && bcQ[c].vertex[X1Y1Z0] == VALUE) { 
			row[l].i = nx-1; row[l].j = ny-1; row[l].k = 0; row[l].c = c; 
			l++;
		}
		if (xs + xm == nx && ys + ym == ny && zs + zm == nz && bcQ[c].vertex[X1Y1Z1] == VALUE && dim ==3) { 
			row[l].i = nx=1; row[l].j = ny-1; row[l].k = nz-1; row[l].c = c; 
			l++;
		}
		
	}
	ierr = MatZeroRowsStencil(K,numBC,row,one,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
		//	ierr = MatZeroRowsStencil(Klhs,numBC,row,zero,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscFree(row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


