/*
 VFFlow_TSFEM.c
 Finite element implementation of  Darcy flow 
 (c) 2010-2013 K. Yoshioka
 */

#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFFlow.h"
#include "VFFlow_TSFEM.h"

#undef __FUNCT__
#define __FUNCT__ "FEMTSFlowSolverInitialize"
extern PetscErrorCode FEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields)
{
	PetscMPIInt    comm_size;
	PetscErrorCode ierr;
	
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&comm_size);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,&ctx->KP);CHKERRQ(ierr);
  ierr = DMCreateMatrix(ctx->daScal,&ctx->KPlhs);CHKERRQ(ierr);		
  ierr = DMCreateMatrix(ctx->daScal,&ctx->JacP);CHKERRQ(ierr);
  ierr = MatZeroEntries(ctx->JacP);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->KPlhs,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = MatSetOption(ctx->JacP,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(ctx->daScal,&ctx->RHSP);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) ctx->RHSP,"RHS of FEM Pressure");CHKERRQ(ierr);	
	ierr = DMCreateGlobalVector(ctx->daScal,&ctx->PFunct);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject)ctx->PFunct,"RHS of FEM TS flow solver");CHKERRQ(ierr);

	ierr = BCPInit(&ctx->bcP[0],ctx);
	ierr = GetFlowProp(&ctx->flowprop,ctx->units,&ctx->resprop,ctx->matprop,ctx,fields);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FEMTSFlowSolverFinalize"
extern PetscErrorCode FEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = MatDestroy(&ctx->KP);CHKERRQ(ierr);
	ierr = MatDestroy(&ctx->KPlhs);CHKERRQ(ierr);	
	ierr = MatDestroy(&ctx->JacP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->RHSP);CHKERRQ(ierr);
	ierr = VecDestroy(&ctx->PFunct);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FlowFEMTSSolve"
extern PetscErrorCode FlowFEMTSSolve(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	PetscInt           temp_step=0;

	PetscFunctionBegin;	
	ierr = DMCreateGlobalVector(ctx->daVFperm,&ctx->Perm);CHKERRQ(ierr);
	ierr = VecSet(ctx->Perm,0.0);CHKERRQ(ierr);
	ierr = VecCopy(fields->vfperm,ctx->Perm);CHKERRQ(ierr);
	
	ierr = TSCreate(PETSC_COMM_WORLD,&ctx->tsP);CHKERRQ(ierr);
	ierr = TSSetDM(ctx->tsP,ctx->daScal);CHKERRQ(ierr);
	ierr = TSSetProblemType(ctx->tsP,TS_LINEAR);CHKERRQ(ierr);
	ierr = TSSetType(ctx->tsP,TSBEULER);CHKERRQ(ierr);	
	
	ierr = TSSetIFunction(ctx->tsP,PETSC_NULL,FormIFunction_P,ctx);CHKERRQ(ierr);
    ierr = TSSetIJacobian(ctx->tsP,ctx->JacP,ctx->JacP,FormIJacobian_P,ctx);CHKERRQ(ierr);
	ierr = TSSetRHSFunction(ctx->tsP,PETSC_NULL,FormFunction_P,ctx);CHKERRQ(ierr);
  
	ierr = TSSetSolution(ctx->tsP,fields->pressure);CHKERRQ(ierr);
	ierr = TSSetInitialTimeStep(ctx->tsP,ctx->current_time,ctx->flowprop.timestepsize);CHKERRQ(ierr);
    ierr = TSSetDuration(ctx->tsP,ctx->maxtimestep,ctx->timevalue);CHKERRQ(ierr);
   
	ierr = TSMonitorSet(ctx->tsP,FEMTSMonitor,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = TSSetFromOptions(ctx->tsP);CHKERRQ(ierr);	

    ierr = TSSolve(ctx->tsP,fields->pressure);CHKERRQ(ierr);
    ierr = TSGetTimeStepNumber(ctx->tsP,&temp_step);CHKERRQ(ierr);
	
	ierr = VecDestroy(&ctx->Perm);CHKERRQ(ierr);
	ierr = TSDestroy(&ctx->tsP);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "FEMTSMonitor"
extern PetscErrorCode FEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec pressure,void* ptr)
{
	PetscErrorCode ierr;
	PetscReal      norm,vmax,vmin;
	MPI_Comm       comm;

	PetscFunctionBegin;
	ierr = VecNorm(pressure,NORM_1,&norm);CHKERRQ(ierr);
	ierr = VecMax(pressure,PETSC_NULL,&vmax);CHKERRQ(ierr);
	ierr = VecMin(pressure,PETSC_NULL,&vmin);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)ts,&comm);CHKERRQ(ierr);
	ierr = PetscPrintf(comm,"timestep %D: time %G, solution norm %G, max %G, min %G\n",timestep,timevalue,norm,vmax,vmin);CHKERRQ(ierr);
		
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormFunction_P"
extern PetscErrorCode FormFunction_P(TS ts,PetscReal t,Vec vec1,Vec Func,void *user)
{
	PetscErrorCode ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = VecCopy(ctx->RHSP,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FormIJacobian_P"
extern PetscErrorCode FormIJacobian_P(TS ts,PetscReal t,Vec pressure,Vec pressuredot,PetscReal shift,Mat Jac,Mat Jacpre,void *user)
{
	PetscErrorCode ierr;
	VFCtx          *ctx=(VFCtx*)user;

	PetscFunctionBegin;
	ierr = MatCopy(ctx->KP,Jacpre,DIFFERENT_NONZERO_PATTERN);
	ierr = MatAXPY(Jacpre,shift,ctx->KPlhs,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Jacpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	if (Jac != Jacpre) {
  	ierr = MatCopy(Jacpre,Jac,DIFFERENT_NONZERO_PATTERN);
		ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}	
	
/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"Jacobian.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = MatView(*Jac,viewer);CHKERRQ(ierr); */

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormTSMatricesnVector_P"
extern PetscErrorCode FormTSMatricesnVector_P(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx)
{
	PetscErrorCode ierr;
	PetscInt       xs,xm,nx;
	PetscInt       ys,ym,ny;
	PetscInt       zs,zm,nz;
	PetscInt       ek,ej,ei;
	PetscInt       i,j,k,l,ii;
	PetscReal      ****perm_array;
	PetscReal      ****coords_array;
	PetscReal      ***RHS_array;
	Vec            RHS_localVec;	
	PetscReal      *RHS_local;
	Vec            perm_local;
	PetscReal      hx,hy,hz;
	PetscReal      *K_local,*Klhs_local;
	PetscInt       nrow = ctx->e3D.nphix*ctx->e3D.nphiy*ctx->e3D.nphiz;
	MatStencil     *row;
	PetscReal      ***source_array;
	Vec            source_local;
	PetscReal      hwx,hwy,hwz;
	
	PetscFunctionBegin;
	ierr = DMDAGetInfo(ctx->daScalCell,PETSC_NULL,&nx,&ny,&nz,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = DMDAGetCorners(ctx->daScalCell,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = MatZeroEntries(K);CHKERRQ(ierr);
	ierr = MatZeroEntries(Klhs);CHKERRQ(ierr);
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

	ierr = PetscMalloc2(nrow*nrow,&K_local,nrow*nrow,&Klhs_local);CHKERRQ(ierr);
	ierr = PetscMalloc2(nrow,&RHS_local,nrow,&row);CHKERRQ(ierr);
	
	for (ek = zs; ek < zs+zm; ek++) {
		for (ej = ys; ej < ys+ym; ej++) {
			for (ei = xs; ei < xs+xm; ei++) {
				hx   = coords_array[ek][ej][ei+1][0]-coords_array[ek][ej][ei][0];
				hy   = coords_array[ek][ej+1][ei][1]-coords_array[ek][ej][ei][1];
				hz   = coords_array[ek+1][ej][ei][2]-coords_array[ek][ej][ei][2];
				ierr = VFCartFEElement3DInit(&ctx->e3D,hx,hy,hz);CHKERRQ(ierr);
				
				ierr = VFFlow_FEM_MatMPAssembly3D_local(Klhs_local,&ctx->flowprop,ek,ej,ei,&ctx->e3D);CHKERRQ(ierr);			
                ierr = VFFlow_FEM_MatKPAssembly3D_local(K_local,&ctx->flowprop,perm_array,ek,ej,ei,&ctx->e3D);					
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							row[l].i = ei+i;row[l].j = ej+j;row[l].k = ek+k;row[l].c = 0;
						}
					}
				}
				ierr = MatSetValuesStencil(K,nrow,row,nrow,row,K_local,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValuesStencil(Klhs,nrow,row,nrow,row,Klhs_local,ADD_VALUES);CHKERRQ(ierr);
				/*Assembling the right hand side vector g*/
				ierr = Flow_Vecg_FEM(RHS_local,&ctx->e3D,ek,ej,ei,ctx->flowprop,perm_array);CHKERRQ(ierr);
				for (l = 0,k = 0; k < ctx->e3D.nphiz; k++) {
					for (j = 0; j < ctx->e3D.nphiy; j++) {
						for (i = 0; i < ctx->e3D.nphix; i++,l++) {
							RHS_array[ek+k][ej+j][ei+i] += RHS_local[l];
						}
					}
				}
				/* Adding source term */
				ierr = VecApplySourceTerms_FEM(RHS_local,source_array,&ctx->e3D,ek,ej,ei,ctx);
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
								ierr = VecApplyWellFlowRate_FEM(RHS_local,ctx->well[ii].Qw,hwx,hwy,hwz);CHKERRQ(ierr);
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
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatApplyPressureBC_FEM(K,Klhs,&ctx->bcP[0]);CHKERRQ(ierr);

    ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Klhs,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		
	ierr = DMDAVecRestoreArray(ctx->daScal,source_local,&source_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&source_local);CHKERRQ(ierr);	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVFperm,perm_local,&perm_array);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daVFperm,&perm_local);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(ctx->daScal,RHS_localVec,&RHS_array);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(ctx->daScal,RHS_localVec,ADD_VALUES,RHS);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(ctx->daScal,&RHS_localVec);CHKERRQ(ierr);

    ierr = VecApplyPressureBC_FEM(RHS,ctx->PresBCArray,&ctx->bcP[0]);CHKERRQ(ierr);
	
	ierr = DMDAVecRestoreArrayDOF(ctx->daVect,ctx->coordinates,&coords_array);CHKERRQ(ierr);

	ierr = PetscFree2(K_local,Klhs_local);CHKERRQ(ierr);
	ierr = PetscFree2(RHS_local,row);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormIFunction_P"
extern PetscErrorCode FormIFunction_P(TS ts,PetscReal t,Vec pressure,Vec pressuredot,Vec Func,void *user)
{
	PetscErrorCode  ierr;
	VFCtx			*ctx=(VFCtx*)user;
	
	PetscFunctionBegin;
	ierr = FormTSMatricesnVector_P(ctx->KP,ctx->KPlhs,ctx->RHSP,ctx);CHKERRQ(ierr);
/*	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"RHSVector.txt",&viewer);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX);CHKERRQ(ierr);
	ierr = VecView(ctx->RHSP,viewer);CHKERRQ(ierr);	*/
	ierr = VecSet(Func,0.0);CHKERRQ(ierr);
	ierr = MatMult(ctx->KP,pressure,Func);CHKERRQ(ierr);	
	ierr = MatMultAdd(ctx->KPlhs,pressuredot,Func,Func);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

