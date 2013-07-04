/*
 VFFlow_SNESFEM.c
 Finite element implementation of  Darcy flow 
 (c) 2010-2013 K. Yoshioka
 */

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_SNESFEM.h"

#undef __FUNCT__
#define __FUNCT__ "FEMSNESFlowSolverInitialize"
extern PetscErrorCode FEMSNESFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->KP);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->KPlhs);CHKERRQ(ierr);		
  ierr = DMCreateMatrix(ctx->daScal,MATAIJ,&ctx->JacP);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->JacP);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	ierr = MatSetOption(ctx->JacP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
	
  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSP);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSPpre);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PFunct);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PrePressure);CHKERRQ(ierr);
	ierr = VecSet(ctx->RHSPpre,0.);CHKERRQ(ierr);	
	ierr = VecSet(ctx->PrePressure,0.);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) ctx->RHSP,"RHS of FEM Pressure");CHKERRQ(ierr);	
	ierr = PetscObjectSetName((PetscObject)ctx->RHSPpre,"RHS of FEM previous pressure");CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->PFunct,"RHS of FEM SNES flow solver");CHKERRQ(ierr);
	
	ierr = SNESCreate(PETSC_COMM_WORLD,&ctx->snesP);CHKERRQ(ierr);
		
	ierr = BCPInit(&ctx->bcP[0],ctx);

	ierr = GetFlowProp(&ctx->flowprop,ctx->units,ctx->resprop);CHKERRQ(ierr);
//	ierr = SETFlowBC(&ctx->bcP[0],&ctx->bcQ[0],ctx->flowcase);CHKERRQ(ierr);	// Currently BCpres is a PetscOption received from command line. Also done in test35
	ierr = ResetSourceTerms(ctx->Source,ctx->flowprop);
//	ierr = SETBoundaryTerms_P(ctx,fields);CHKERRQ(ierr); // Set fields.pressure according to BC & IC
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FEMSNESFlowSolverFinalize"
extern PetscErrorCode FEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
    ierr = MatDestroy(&ctx->KP);CHKERRQ(ierr);
    ierr = MatDestroy(&ctx->KPlhs);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->JacP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PrePressure);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSPpre);CHKERRQ(ierr);	
	ierr = VecDestroy(&ctx->PFunct);CHKERRQ(ierr);
	ierr = SNESDestroy(&ctx->snesP);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowFEMSNESSolve"
extern PetscErrorCode FlowFEMSNESSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;

	PetscFunctionBegin;	
	ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Perm);CHKERRQ(ierr);
	ierr = VecSet(ctx->Perm,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->vfperm,ctx->Perm);CHKERRQ(ierr);
  ierr = VecCopy(fields->pressure,ctx->PrePressure);CHKERRQ(ierr);
	
	ierr = SNESSetFunction(ctx->snesP,ctx->PFunct,FormSNESIFunction_P,ctx);CHKERRQ(ierr);
  ierr = SNESSetJacobian(ctx->snesP,ctx->JacP,ctx->JacP,FormSNESIJacobian_P,ctx);CHKERRQ(ierr);
	ierr = SNESMonitorSet(ctx->snesP,FEMSNESMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = SNESSolve(ctx->snesP,PETSC_NULL,fields->pressure);CHKERRQ(ierr);
  ierr = VecCopy(ctx->RHSP,ctx->RHSPpre);CHKERRQ(ierr);	
		
	ierr = VecDestroy(&ctx->Perm);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESIJacobian_P"
extern PetscErrorCode FormSNESIJacobian_P(SNES snes,Vec pressure,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user)
{
	PetscErrorCode ierr;
	VFCtx          *ctx=(VFCtx*)user;

	PetscFunctionBegin;
	*str = DIFFERENT_NONZERO_PATTERN;
	ierr = MatZeroEntries(*Jac);CHKERRQ(ierr);
	ierr = MatCopy(ctx->KP,*Jac,*str);
	ierr = MatAssemblyBegin(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (*Jac != *Jacpre) {
		ierr = MatAssemblyBegin(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(*Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}	

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESMatricesnVector_P"
extern PetscErrorCode FormSNESMatricesnVector_P(Mat Kneu,Mat Kalt,Vec RHS,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm;
	PetscInt       ys,ym;
	PetscInt       zs,zm;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l,ii;
	PetscReal      ****perm_array;
	PetscReal      ****coords_array;
	PetscReal      ***RHS_array;
	Vec            RHS_localVec;	
	PetscReal      *RHS_local;
	Vec            perm_local;
	PetscReal      hx,hy,hz;
	PetscReal      *K_local,*M_local;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row;
	PetscReal      ***source_array;
	Vec            source_local;
	Mat            M,K;
	PetscReal	     theta,timestepsize;	
	PetscReal      time_theta,time_one_minus_theta;	
	PetscReal      hwx,hwy,hwz;
	MatStructure   flg;
	
	PetscFunctionBegin;
	flg = SAME_NONZERO_PATTERN;
	theta = ctx->flowprop.theta;
	timestepsize = ctx->flowprop.timestepsize;
	time_theta = theta * timestepsize;
	time_one_minus_theta = -1.*(1.-theta)*timestepsize;
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = MatZeroEntries(Kneu);CHKERRQ(ierr);
	ierr = MatZeroEntries(Kalt);CHKERRQ(ierr);
	ierr = MatDuplicate(Kneu,MAT_COPY_VALUES,&M);CHKERRQ(ierr);
	ierr = MatDuplicate(Kneu,MAT_COPY_VALUES,&K);CHKERRQ(ierr);	
	ierr = VecSet(RHS,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
	ierr = VecSet(RHS_localVec,0.);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daScal,ctx->Source,INSERT_VALUES,source_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMGetLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(ctx->daVFperm,ctx->Perm,INSERT_VALUES,perm_local);CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);	
	ierr = PetscMalloc2(nrow*nrow,PetscReal,&K_local,nrow*nrow,PetscReal,&M_local);CHKERRQ(ierr);
	ierr = PetscMalloc2(nrow,PetscReal,&RHS_local,nrow,MatStencil,&row);CHKERRQ(ierr);
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = CartFE_Element3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);

				ierr = VFFlow_FEM_MatMPAssembly3D_local(M_local,&ctx->flowprop,ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);
                ierr = VFFlow_FEM_MatKPAssembly3D_local(K_local,&ctx->flowprop,perm_array,ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);			
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 0;
						}
					}
				}	
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(M,nrow,row,nrow,row,M_local,ADD_VALUES);CHKERRQ(ierr);				
				/*Assembling the right hand side vector g*/
				ierr = Flow_Vecg(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
						}
					}
				}
				/* Adding source term */
				ierr = VecApplySourceTerms(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
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
								ierr = VecApplyWellFlowRate(RHS_local,ctx->well[ii].Qw,hwx,hwy,hwz);CHKERRQ(ierr);
								for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
									for (j = 0; j < ctx->e3D.nphiy; j++) {
										for (i = 0; i < ctx->e3D.nphix; i++,l++) {
											if(ctx->well[ii].type == INJECTOR){
												RHS_array[ek+k][ej+j][ei+i] -= RHS_local[l];
											}
											else if(ctx->well[ii].type == PRODUCER)
											{
												RHS_array[ek+k][ej+j][ei+i] -= -1*RHS_local[l];
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
	
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 

	ierr = MatAXPY(Kneu,time_theta,K,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Kneu,1.0,M,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Kalt,time_one_minus_theta,K,flg);CHKERRQ(ierr);
	ierr = MatAXPY(Kalt,1.0,M,flg);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Kalt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Kalt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Kneu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Kneu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatApplyPressureBC_FEM(Kneu,Kalt,&ctx->bcP[0]);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Kneu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Kneu,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Kalt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Kalt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

	ierr = PetscFree2(K_local,M_local);CHKERRQ(ierr);
	ierr = PetscFree2(RHS_local,row);CHKERRQ(ierr);	
	ierr = MatDestroy(&M);CHKERRQ(ierr);
	ierr = MatDestroy(&K);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormSNESIFunction_P"
extern PetscErrorCode FormSNESIFunction_P(SNES snes,Vec pressure,Vec Func,void *user)
{
	PetscErrorCode  ierr;
	VFCtx			     *ctx=(VFCtx*)user;
	Vec             VecRHS;
	PetscReal	      theta,timestepsize;
	PetscReal		    dt_dot_theta,dt_dot_one_minus_theta;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	timestepsize = ctx->flowprop.timestepsize;
	theta = ctx->flowprop.theta;
	dt_dot_theta = timestepsize*theta;
	dt_dot_one_minus_theta = timestepsize*(1.-theta);
	ierr = VecDuplicate(ctx->RHSP,&VecRHS);CHKERRQ(ierr);
	ierr = FormSNESMatricesnVector_P(ctx->KP,ctx->KPlhs,ctx->RHSP,ctx);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSP,VecRHS);CHKERRQ(ierr);
	ierr = VecAXPBY(VecRHS,dt_dot_one_minus_theta,dt_dot_theta,ctx->RHSPpre);CHKERRQ(ierr);	
	// VecCopy RHSP RHSPpre?
	ierr = MatMultAdd(ctx->KPlhs,pressure,VecRHS,VecRHS);CHKERRQ(ierr);	
	ierr = VecApplyPressureBC_FEM(VecRHS,ctx->PresBCArray,&ctx->bcP[0]);CHKERRQ(ierr);

/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHS.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(VecRHS,viewer);CHKERRQ(ierr);	
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Mat.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = MatView(ctx->KP,viewer);CHKERRQ(ierr);	*/
	
	ierr = MatMult(ctx->KP,pressure,Func);CHKERRQ(ierr);	
	ierr = VecAXPY(Func,-1.0,VecRHS);CHKERRQ(ierr);
	ierr = VecDestroy(&VecRHS);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}	




